// Enable this to output polygon debug information
//#define WITH_DEBUG
using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Unity.Mathematics;
using UnityEngine;

namespace AggroBird.StraightSkeleton
{
    public sealed class StraightSkeletonException : Exception
    {
        public StraightSkeletonException(string msg) : base(msg)
        {

        }
    }

    public sealed class StraightSkeleton : IReadOnlyList<IReadOnlyList<float2>>
    {
        private readonly List<List<float2>> polygons = new List<List<float2>>();

        public int Count { get; private set; }
        public IReadOnlyList<float2> this[int index] => polygons[index];

        public IEnumerator<IReadOnlyList<float2>> GetEnumerator()
        {
            for (int i = 0; i < Count; i++)
            {
                yield return polygons[i];
            }
        }
        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public void Clear()
        {
            Count = 0;

#if WITH_DEBUG
            debugOutput.Clear();
#endif
        }

        internal List<float2> GetBuffer()
        {
            List<float2> buffer;
            if (polygons.Count == Count)
            {
                buffer = new List<float2>();
                polygons.Add(buffer);
                Count++;
            }
            else
            {
                buffer = polygons[Count++];
                buffer.Clear();
            }
            return buffer;
        }

#if WITH_DEBUG
        private readonly List<(float2, float2, Color)> debugOutput = new List<(float2, float2, Color)>();
        public IReadOnlyList<(float2 from, float2 to, Color color)> DebugOutput => debugOutput;
        internal void AddDebugLine(float2 from, float2 to, Color color)
        {
            debugOutput.Add((from, to, color));
        }
#else
        public IReadOnlyList<(float2 from, float2 to, Color color)> DebugOutput => Array.Empty<(float2, float2, Color)>();
#endif
    }

    public sealed class StraightSkeletonGenerator
    {
        private const float HighPrecisionEpsilon = math.EPSILON;
        private const float MediumPrecisionEpsilon = 0.00001f;
        private const float LowPrecisionEpsilon = 0.001f;

#if WITH_DEBUG
        private static readonly Color[] DebugColors = new Color[]
        {
            Color.red,
            Color.green,
            Color.blue,
            Color.yellow,
            Color.magenta,
            Color.cyan,
            new Color32(255, 174, 201, 255), // Pink
            new Color32(255, 127, 39, 255), // Orange
            new Color32(163, 73, 164, 255), // Purple
        };
#endif


        private enum PolygonSide
        {
            Unknown,
            Left = 1,
            Right = 2,
        }

        private struct PolygonVertex
        {
            public PolygonVertex(float2 position, int pass = 0)
            {
                index = -1;

                this.position = position;
                this.pass = pass;

                prevPolyVert = -1;
                nextPolyVert = -1;
            }

            public int index;

            // The index of the pass this vertex was inserted
            public int pass;

            // 2D position of the vertex
            public float2 position;

            // Neighbours in the polygon linked list
            public int prevPolyVert;
            public int nextPolyVert;


            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void Link(ref PolygonVertex prev, ref PolygonVertex next)
            {
                LinkPrev(ref prev);
                LinkNext(ref next);
            }
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void LinkPrev(ref PolygonVertex prev)
            {
                prev.nextPolyVert = index;
                prevPolyVert = prev.index;
            }
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void LinkNext(ref PolygonVertex next)
            {
                next.prevPolyVert = index;
                nextPolyVert = next.index;
            }
        }

        private struct ChainVertex
        {
            public int index;

            // The index of the pass this vertex was inserted
            public int pass;

            // 2D position of the vertex
            public float2 position;

            // Direction of the segment to the next vertex (normalized)
            public float2 direction;
            // Distance to the next vertex
            public float length;

            // Angle of the vertex between next and previous
            public float angle;
            // Direction of the bisector (normalized)
            public float2 bisector;
            // Speed at which the bisector shrinks
            public float velocity;

            public bool IsReflex { get; private set; }

            // Neighbours in the chain
            public int prevChainVert;
            public int nextChainVert;

            // Neighbouring polygons on this axis
            public int lhsPolyVert;
            public int rhsPolyVert;


            // Recalculate direction and length
            public void RecalculateSegment(float2 nextPosition)
            {
                direction = nextPosition - position;
                length = math.length(direction);
                if (length <= HighPrecisionEpsilon)
                {
                    direction.x = direction.y = 0;
                    length = 0;
                }
                else
                {
                    direction /= length;
                }
            }
            // Recalculate bisector and angle
            public void RecalculateBisector(float2 prevDirection)
            {
                float2 a = direction;
                float2 b = -prevDirection;

                float signedAngle = math.atan2(b.x * a.y - b.y * a.x, math.dot(b, a));
                float PI2 = math.PI * 2;
                angle = (signedAngle + PI2) % PI2;
                float angleHalf = angle * 0.5f;
                bisector = Rotate(a, angleHalf);
                velocity = math.sin(angleHalf);
                IsReflex = angle > math.PI;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void Link(ref ChainVertex prev, ref ChainVertex next)
            {
                prev.nextChainVert = index;
                prevChainVert = prev.index;
                next.prevChainVert = index;
                nextChainVert = next.index;
            }


            public override string ToString()
            {
                return $"chain vertex #{index}";
            }
        }

        // Chain of vertices.
        // Vertices are not guaranteed to be in order, so they are stored as a linked list.
        /*private class Chain : IEnumerable<ChainVertex>
        {

            public ChainVertex First => vertices[0];
            public int Count => vertices.Count;
            public ChainVertex this[int idx] => vertices[idx];
            private readonly List<ChainVertex> vertices = new List<ChainVertex>();

            public void Add(ChainVertex vertex)
            {
                vertices.Add(vertex);
            }
            public void Clear()
            {
                vertices.Clear();
            }

            // Enumerators for the linked lists
            public IEnumerator<ChainVertex> GetEnumerator()
            {
                if (vertices.Count == 0) yield break;
                var start = First;
                var current = start;
                do
                {
                    yield return current;
                    current = current.nextChainVert;
                }
                while (!current.Equals(start));
            }
            IEnumerator IEnumerable.GetEnumerator()
            {
                return GetEnumerator();
            }



            public override string ToString()
            {
                if (vertices.Count == 0)
                {
                    return "empty chain";
                }
                else
                {
                    return $"chain of {First} ({vertices.Count}) vertices";
                }
            }
        }*/

        private static void EnsureCapacity<T>(ref T[] buffer, int size, bool copyData)
        {
            int currentCapacity = buffer.Length;
            if (size > currentCapacity)
            {
                int newCapacity = currentCapacity;
            IncreaseCapacity:
                newCapacity <<= 1;
                if (newCapacity <= 0)
                {
                    throw new StraightSkeletonException("Buffer allocation caused an overflow");
                }
                if (size > newCapacity)
                {
                    goto IncreaseCapacity;
                }

                T[] newBuffer = new T[newCapacity];
                if (copyData)
                {
                    Array.Copy(buffer, newBuffer, currentCapacity);
                }
                buffer = newBuffer;
            }
        }

        private ChainVertex[] chainVertices = new ChainVertex[16];
        private int chainVertexCount = 0;
        private ref ChainVertex AddChainVertex(ChainVertex chainVertex)
        {
            chainVertex.index = chainVertexCount;
            EnsureCapacity(ref chainVertices, chainVertexCount + 1, true);
            ref ChainVertex result = ref chainVertices[chainVertexCount++];
            result = chainVertex;
            return ref result;
        }

        private readonly List<int> activeChains = new List<int>();
        private readonly List<int> newChains = new List<int>();

        private PolygonVertex[] polygonVertices = new PolygonVertex[32];
        private int polygonVertexCount = 0;
        private ref PolygonVertex AddPolygonVertex(PolygonVertex polygonVertex)
        {
            polygonVertex.index = polygonVertexCount;
            EnsureCapacity(ref polygonVertices, polygonVertexCount + 1, true);
            ref PolygonVertex result = ref polygonVertices[polygonVertexCount++];
            result = polygonVertex;
            return ref result;
        }

        private readonly List<int> incidentVertices = new List<int>();
        private readonly List<int> unresolvedChains = new List<int>();
        private readonly List<int> resolvedChains = new List<int>();


        // Allow reusage of straight skeleton buffer data
        public void Generate(IReadOnlyList<float2> points, StraightSkeleton output)
        {
            output.Clear();

            if (points == null || points.Count < 3)
            {
                throw new StraightSkeletonException("Invalid input polygon provided");
            }

            int pointCount = points.Count;

            // Create initial chain
            chainVertexCount = pointCount;
            EnsureCapacity(ref chainVertices, chainVertexCount, false);
            for (int i = 0; i < chainVertexCount; i++)
            {
                chainVertices[i] = new ChainVertex { index = i, position = points[i] };
            }
            for (int i = 0; i < chainVertexCount; i++)
            {
                int j = (i + 1) % chainVertexCount;
                ref ChainVertex prev = ref chainVertices[i];
                ref ChainVertex next = ref chainVertices[j];
                prev.nextChainVert = j;
                next.prevChainVert = i;
                prev.RecalculateSegment(next.position);
            }
            for (int i = 0; i < chainVertexCount; i++)
            {
                int j = (i + 1) % chainVertexCount;
                ref ChainVertex prev = ref chainVertices[i];
                ref ChainVertex next = ref chainVertices[j];
                next.RecalculateBisector(prev.direction);
            }
            activeChains.Clear();
            activeChains.Add(0);

            // Create initial starting polygons
            polygonVertexCount = chainVertexCount * 2;
            EnsureCapacity(ref polygonVertices, polygonVertexCount, false);
            for (int i = 0, idx = 0; i < chainVertexCount; i++)
            {
                int j = (i + 1) % chainVertexCount;
                ref ChainVertex next = ref chainVertices[j];
                ref ChainVertex prev = ref chainVertices[i];
                int lhsIdx = idx++, rhsIdx = idx++;
                PolygonVertex lhs = new PolygonVertex(next.position);
                PolygonVertex rhs = new PolygonVertex(prev.position);
                lhs.index = lhsIdx;
                rhs.index = rhsIdx;
                lhs.prevPolyVert = rhsIdx;
                rhs.nextPolyVert = lhsIdx;
                next.lhsPolyVert = lhsIdx;
                prev.rhsPolyVert = rhsIdx;
                polygonVertices[lhsIdx] = lhs;
                polygonVertices[rhsIdx] = rhs;
            }

#if WITH_DEBUG
            for (int i = 0; i < chainVertexCount; i++)
            {
                int j = (i + 1) % chainVertexCount;
                output.AddDebugLine(chainVertices[i].position, chainVertices[j].position, Color.white);
                output.AddDebugLine(chainVertices[i].position, chainVertices[i].position + chainVertices[i].bisector * 0.1f, Color.blue);
            }
#endif

            int pass = 0;
            while (activeChains.Count > 0)
            {
                if (pass++ > pointCount)
                {
                    throw new StraightSkeletonException("Max polygon shrink iteration count reached");
                }

                newChains.Clear();
                for (int i = 0; i < activeChains.Count; i++)
                {
                    int activeChain = activeChains[i];

                    if (!CalculateShortestShrinkDistance(activeChain, out float distance))
                    {
                        throw new StraightSkeletonException("Failed to find distance shortest shrink distance");
                    }

                    ApplyShrinkDistance(activeChain, distance);

                    ProcessIntersectionEvents(activeChain, pass);
                }

                activeChains.Clear();
                activeChains.AddRange(newChains);

#if WITH_DEBUG
                foreach (var chain in activeChains)
                {
                    int vertIdx = chain;
                    do
                    {
                        ref ChainVertex vert = ref chainVertices[vertIdx];
                        ref ChainVertex next = ref chainVertices[vert.nextChainVert];
                        output.AddDebugLine(vert.position, next.position, Color.grey);
                        output.AddDebugLine(vert.position, vert.position + vert.bisector * 0.1f, Color.blue);
                        vertIdx = vert.nextChainVert;
                    }
                    while (vertIdx != chain);
                }
#endif
            }

#if WITH_DEBUG
            List<float2> verts = new List<float2>();
            List<float2> bisectors = new List<float2>();
            bool isValid = true;
            for (int i = 0; i < pointCount; i++)
            {
                verts.Clear();
                bisectors.Clear();

                int first = i * 2;
                int index = first;
                int iter = 0;
                bool hasInfiniteLoop = false;

                // Try and find the initial vertex (in case of unconnected loops)
                do
                {
                    if (iter++ > pointCount * 100)
                    {
                        Debug.LogError($"polygon #{i} contains an infinite loop");
                        hasInfiniteLoop = true;
                        isValid = false;
                        break;
                    }
                    ref PolygonVertex vertex = ref polygonVertices[index];
                    if (vertex.prevPolyVert == -1)
                    {
                        break;
                    }
                    index = vertex.prevPolyVert;
                }
                while (!index.Equals(first));
                first = index;

                // Fetch all vertices
                iter = 0;
                do
                {
                    if (iter++ > pointCount * 100)
                    {
                        if (!hasInfiniteLoop)
                        {
                            Debug.LogError($"polygon #{i} contains an infinite loop");
                            hasInfiniteLoop = true;
                        }
                        isValid = false;
                        break;
                    }
                    ref PolygonVertex vertex = ref polygonVertices[index];
                    verts.Add(vertex.position);
                    if (vertex.nextPolyVert == -1)
                    {
                        Debug.LogError($"polygon #{i} contains null connections");
                        isValid = false;
                        break;
                    }
                    index = vertex.nextPolyVert;
                }
                while (index != first);

                if (verts.Count >= 3)
                {
                    // Generate bisectors for the debug shapes
                    for (int j = 0; j < verts.Count; j++)
                    {
                        int k = (j + 1) % verts.Count;
                        int l = (j + verts.Count - 1) % verts.Count;
                        float2 a = math.normalize(verts[k] - verts[j]);
                        float2 b = -math.normalize(verts[j] - verts[l]);
                        float signedAngle = math.atan2(b.x * a.y - b.y * a.x, math.dot(b, a));
                        float PI2 = math.PI * 2;
                        float angle = (signedAngle + PI2) % PI2;
                        float angleHalf = angle * 0.5f;
                        float2 bisector = Rotate(a, angleHalf);
                        float velocity = math.sin(angleHalf);
                        bisectors.Add(bisector * (0.01f / velocity));
                    }

                    // Output polygons
                    Color debugColor = DebugColors[i % DebugColors.Length];
                    for (int j = 0; j < verts.Count; j++)
                    {
                        int k = (j + 1) % verts.Count;
                        output.AddDebugLine(verts[j] + bisectors[j], verts[k] + bisectors[k], debugColor);
                    }
                }
                else
                {
                    Debug.LogError($"Polygon {i} does not contain enough vertices");
                    isValid = false;
                }
            }
            if (!isValid) return;
#endif

            // Output result
            for (int i = 0; i < pointCount; i++)
            {
                List<float2> buffer = output.GetBuffer();
                int first = i * 2;
                int index = first;
                do
                {
                    ref PolygonVertex vertex = ref polygonVertices[index];
                    buffer.Add(vertex.position);
                    index = vertex.nextPolyVert;
                }
                while (index != first);
            }
        }
        public StraightSkeleton Generate(IReadOnlyList<float2> points)
        {
            StraightSkeleton output = new StraightSkeleton();
            Generate(points, output);
            return output;
        }

        // Calculates the shortest shrink distance at which point a collision will occur
        private bool CalculateShortestShrinkDistance(int chain, out float distance)
        {
            distance = float.PositiveInfinity;

            int vertIdx = chain;
            do
            {
                ref ChainVertex vert = ref chainVertices[vertIdx];
                ref ChainVertex next = ref chainVertices[vert.nextChainVert];

                // Project the two bisectors onto each other and find the shortest collapse distance (edge event)
                float2 n0 = new float2(-vert.bisector.y, vert.bisector.x);
                float2 n1 = new float2(next.bisector.y, -next.bisector.x);

                float entry0 = math.dot(next.position - vert.position, n1) / math.dot(vert.bisector, n1);
                float entry1 = math.dot(vert.position - next.position, n0) / math.dot(next.bisector, n0);

                float entry = math.min(entry0 * vert.velocity, entry1 * next.velocity);
                if (!float.IsInfinity(entry) && entry > 0)
                {
                    if (entry < distance)
                    {
                        distance = entry;
                    }
                }

                if (vert.IsReflex)
                {
                    // Find the earliest segment intersect with a reflex vertex (split event)
                    int segIdx = vertIdx;
                    do
                    {
                        ref ChainVertex segBeg = ref chainVertices[segIdx];
                        ref ChainVertex segEnd = ref chainVertices[segBeg.nextChainVert];
                        if (segIdx != vertIdx && segBeg.nextChainVert != vertIdx)
                        {
                            float2 perp = new float2(segBeg.direction.y, -segBeg.direction.x);
                            if (TracePlane(segBeg.position, perp, vert.position, vert.bisector, out entry))
                            {
                                // Calculate time of collision (magic)
                                float sin = math.sin(AngleBetween(vert.bisector, segBeg.direction));
                                entry /= (1 / sin + 1 / vert.velocity);
                                if (!float.IsInfinity(entry) && entry > 0)
                                {
                                    // Grow the segment along the two bisectors of its vertices and project
                                    // the collision point to make sure this vertex will actually collide with
                                    // this segment in the future.
                                    float2 s0 = segBeg.position + segBeg.bisector * (entry / segBeg.velocity);
                                    float2 s1 = segEnd.position + segEnd.bisector * (entry / segEnd.velocity);
                                    float length = math.dot(s1 - s0, segBeg.direction);
                                    if (length > 0)
                                    {
                                        float2 point = vert.position + vert.bisector * (entry / vert.velocity);
                                        float project = math.dot(point - s0, segBeg.direction);
                                        if (project >= -MediumPrecisionEpsilon && project <= length + MediumPrecisionEpsilon)
                                        {
                                            if (entry < distance)
                                            {
                                                distance = entry;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        segIdx = segBeg.nextChainVert;
                    }
                    while (segIdx != vertIdx);
                }

                vertIdx = vert.nextChainVert;
            }
            while (vertIdx != chain);

            return distance != float.PositiveInfinity;
        }

        // Shrink the polygon
        private void ApplyShrinkDistance(int chain, float distance)
        {
            int vertIdx = chain;
            do
            {
                ref ChainVertex vert = ref chainVertices[vertIdx];
                vert.position += vert.bisector * (distance / vert.velocity);
                vertIdx = vert.nextChainVert;
            }
            while (vertIdx != chain);
            do
            {
                ref ChainVertex prev = ref chainVertices[vertIdx];
                ref ChainVertex next = ref chainVertices[prev.nextChainVert];
                prev.RecalculateSegment(next.position);
                vertIdx = prev.nextChainVert;
            }
            while (vertIdx != chain);
        }

        // Check for any intersections (results in one or more new chains)
        private void ProcessIntersectionEvents(int chain, int pass)
        {
            // This phase scans over all the vertices and checks for overlapping (incident) vertices.
            // In the case of a split event, it will emit a new vertex parallel to the segment so that
            // it can be picked up later by the split stage.
            int vertIdx = chain;
            do
            {
                ref ChainVertex vert = ref chainVertices[vertIdx];

                // No need to resolve intersections with inserted split polygons
                if (vert.pass != pass)
                {
                    // Find closest intersect
                    IntersectResult bestResult = IntersectResult.None;
                    int segIdx = vertIdx;
                    do
                    {
                        ref ChainVertex seg = ref chainVertices[segIdx];
                        if (segIdx != vertIdx && seg.nextChainVert != vertIdx)
                        {
                            AppendIntersect(ref vert, ref seg, ref chainVertices[seg.nextChainVert], LowPrecisionEpsilon, ref bestResult);
                        }

                        segIdx = seg.nextChainVert;
                    }
                    while (segIdx != vertIdx);

                    if (bestResult)
                    {
                        // Snap vertex to intersect location
                        vert.position = bestResult.point;
                        RecalculateSegments(ref vert);

                        if (!bestResult.IsIncident)
                        {
                            // Insert new vertex on segment
                            ref ChainVertex segBeg = ref chainVertices[bestResult.segment];
                            ref ChainVertex segEnd = ref chainVertices[segBeg.nextChainVert];
                            ref ChainVertex insert = ref AddChainVertex(new ChainVertex { position = bestResult.point, pass = pass });
                            insert.Link(ref segBeg, ref segEnd);
                            RecalculateSegments(ref insert);
                            // Create an empty polygon vertex at the split position
                            insert.lhsPolyVert = insert.rhsPolyVert = AddPolygonVertex(new PolygonVertex(bestResult.point, pass)).index;
                        }
                    }
                }

                vertIdx = vert.nextChainVert;
            }
            while (vertIdx != chain);

            // This stage takes all incident vertices and splits them up into seperate chains.
            // It also emits polygon vertices for the skeleton in the case of an edge or split event.
            resolvedChains.Clear();
            unresolvedChains.Clear();
            unresolvedChains.Add(chain);
        TraceNext:
            while (unresolvedChains.Count > 0)
            {
                int lastIndex = unresolvedChains.Count - 1;
                int unresolved = unresolvedChains[lastIndex];
                unresolvedChains.RemoveAt(lastIndex);

                vertIdx = unresolved;
                do
                {
                    ref ChainVertex vert = ref chainVertices[vertIdx];

                    // No need to resolve vertices that we have already checked for incident in this pass
                    if (vert.pass != pass)
                    {
                        vert.pass = pass;

                        incidentVertices.Clear();

                        // Collect all vertices that are incident to the current vertex
                        int otherIdx = vertIdx;
                        do
                        {
                            ref ChainVertex other = ref chainVertices[otherIdx];
                            if (vertIdx != otherIdx && IncidentVertices(vert.position, other.position, LowPrecisionEpsilon * 2))
                            {
                                other.pass = pass;

                                incidentVertices.Add(otherIdx);
                            }

                            otherIdx = other.nextChainVert;
                        }
                        while (otherIdx != vertIdx);

                        if (incidentVertices.Count > 0)
                        {
                            // Resolve the current vertex as well
                            incidentVertices.Add(vertIdx);

                            // Split into subshapes
                            int incidentCount = incidentVertices.Count;
                            for (int i = 0; i < incidentCount; i++)
                            {
                                ref ChainVertex prevIncident = ref chainVertices[incidentVertices[i]];
                                ref ChainVertex nextIncident = ref chainVertices[incidentVertices[(i + 1) % incidentCount]];

                                // If there are vertices between the previous and next incident vertices,
                                // we have to split them off into a new chain
                                if (prevIncident.nextChainVert != nextIncident.index)
                                {
                                    ref ChainVertex insert = ref AddChainVertex(new ChainVertex { position = vert.position, pass = pass });

                                    int nextPolygon = prevIncident.rhsPolyVert;
                                    if (polygonVertices[nextPolygon].pass != pass)
                                    {
                                        // Insert a polygon vertex on the right side of the left polygon
                                        ref PolygonVertex lhs = ref AddPolygonVertex(new PolygonVertex(vert.position, pass));
                                        lhs.LinkNext(ref polygonVertices[nextPolygon]);
                                        insert.rhsPolyVert = lhs.index;
                                    }
                                    else
                                    {
                                        // If the polygon point inserted through a split event, we dont have to insert a new vertex
                                        insert.rhsPolyVert = nextPolygon;
                                    }

                                    int prevPolygon = nextIncident.lhsPolyVert;
                                    if (polygonVertices[prevPolygon].pass != pass)
                                    {
                                        // Insert a polygon vertex on the left side of the right polygon
                                        ref PolygonVertex rhs = ref AddPolygonVertex(new PolygonVertex(vert.position, pass));
                                        rhs.LinkPrev(ref polygonVertices[prevPolygon]);
                                        insert.lhsPolyVert = rhs.index;
                                    }
                                    else
                                    {
                                        // If the polygon point inserted through a split event, we dont have to insert a new vertex
                                        insert.lhsPolyVert = prevPolygon;
                                    }

                                    insert.Link(ref chainVertices[nextIncident.prevChainVert], ref chainVertices[prevIncident.nextChainVert]);
                                    RecalculateSegments(ref insert);
                                    unresolvedChains.Add(insert.index);
                                }
                                else
                                {
                                    // Close into a triangle
                                    ref PolygonVertex end = ref AddPolygonVertex(new PolygonVertex(vert.position, pass));
                                    end.Link(ref polygonVertices[nextIncident.lhsPolyVert], ref polygonVertices[prevIncident.rhsPolyVert]);
                                }
                            }

                            goto TraceNext;
                        }
                    }

                    vertIdx = vert.nextChainVert;
                }
                while (vertIdx != unresolved);

                resolvedChains.Add(unresolved);
            }

            // Final stage: rebuild new chains from the split results
            for (int i = 0; i < resolvedChains.Count; i++)
            {
                int vertCount = 0;
                int first = resolvedChains[i];
                vertIdx = first;
                do
                {
                    ref ChainVertex vert = ref chainVertices[vertIdx];
                    vertIdx = vert.nextChainVert;
                    vertCount++;
                }
                while (vertIdx != first);

                if (vertCount > 2)
                {
                    // New valid chain
                    newChains.Add(first);
                }
                else if (vertCount == 2)
                {
                    // If we have only two vertices we can collapse the segment and close the polygons
                    ref ChainVertex prev = ref chainVertices[first];
                    ref ChainVertex next = ref chainVertices[prev.nextChainVert];
                    ref PolygonVertex prevLhs = ref polygonVertices[prev.lhsPolyVert];
                    ref PolygonVertex nextLhs = ref polygonVertices[next.lhsPolyVert];
                    ref PolygonVertex prevRhs = ref polygonVertices[prev.rhsPolyVert];
                    ref PolygonVertex nextRhs = ref polygonVertices[next.rhsPolyVert];
                    if (prevLhs.nextPolyVert == -1 && nextRhs.prevPolyVert == -1)
                    {
                        prevLhs.LinkNext(ref nextRhs);
                    }
                    if (prevRhs.prevPolyVert == -1 && nextLhs.nextPolyVert == -1)
                    {
                        prevRhs.LinkPrev(ref nextLhs);
                    }
                }
            }
        }


        // Recalculate a node and its neighbours
        private void RecalculateSegments(ref ChainVertex vert)
        {
            ref ChainVertex prev = ref chainVertices[vert.prevChainVert];
            ref ChainVertex next = ref chainVertices[vert.nextChainVert];

            prev.RecalculateSegment(vert.position);
            vert.RecalculateSegment(next.position);

            prev.RecalculateBisector(chainVertices[prev.prevChainVert].direction);
            vert.RecalculateBisector(prev.direction);
            next.RecalculateBisector(vert.direction);
        }


        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float AngleBetween(float2 from, float2 to)
        {
            float sqrt = math.sqrt(math.lengthsq(from) * math.lengthsq(to));
            if (sqrt < 1e-15f) return 0f;
            return math.acos(math.clamp(math.dot(from, to) / sqrt, -1f, 1f));
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float2 Rotate(float2 vec, float angle)
        {
            float xr = math.sin(-angle);
            float yr = math.cos(-angle);
            float x = vec.x * yr - vec.y * xr;
            float y = vec.x * xr + vec.y * yr;
            return new float2(x, y);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static bool IncidentVertices(float2 a, float2 b, float accuracy)
        {
            return math.lengthsq(a - b) < accuracy * accuracy;
        }

        // Trace a plane (position, normal) from a ray (origin, direction).
        // Returns false if the ray starts below the plane or the ray is facing away from the plane.
        private static bool TracePlane(float2 position, float2 normal, float2 origin, float2 direction, out float entry)
        {
            if (math.dot(origin - position, normal) >= 0)
            {
                float denominator = math.dot(direction, normal);
                if (math.abs(denominator) >= HighPrecisionEpsilon)
                {
                    entry = math.dot(position - origin, normal) / denominator;
                    return entry >= 0;
                }
            }
            entry = 0;
            return false;
        }

        private struct IntersectResult
        {
            public static readonly IntersectResult None = new IntersectResult { distanceSqr = float.PositiveInfinity, segment = -1 };

            public float2 point;
            public float distanceSqr;
            public int segment;

            public bool IsIncident => segment == -1;

            public static implicit operator bool(IntersectResult intersectResult)
            {
                return intersectResult.distanceSqr != float.PositiveInfinity;
            }
        }

        // Calculate the intersection of a vertex on a segment, and if the distance is closer than the previous, update the result.
        // This function prefers incident results, to reduce the amount of tiny segment intersections.
        private void AppendIntersect(ref ChainVertex vert, ref ChainVertex segBeg, ref ChainVertex segEnd, float accuracy, ref IntersectResult result)
        {
            float accuracySqr = accuracy * accuracy;

            // See if we are anywhere near the segment origin
            float dist = math.distancesq(vert.position, segBeg.position);
            if (dist <= accuracySqr && dist < result.distanceSqr)
            {
                result.point = segBeg.position;
                result.distanceSqr = dist;
                result.segment = -1;
            }

            // See if we are anywhere near the segment end
            dist = math.distancesq(vert.position, segEnd.position);
            if (dist <= accuracySqr && dist < result.distanceSqr)
            {
                result.point = segEnd.position;
                result.distanceSqr = dist;
                result.segment = -1;
            }

            // Project the point onto the segment (if we are currently not incident)
            if (segBeg.length > 0 && (!result || !result.IsIncident))
            {
                // Project on the segment along the length
                float2 relative = vert.position - segBeg.position;
                float project = math.dot(relative, segBeg.direction);
                if (project >= 0 && project <= segBeg.length)
                {
                    // Project on the segment along the width (with an error margin)
                    float2 perp = new float2(segBeg.direction.y, -segBeg.direction.x);
                    float dot = math.dot(perp, relative);
                    if (math.abs(dot) <= accuracy)
                    {
                        // The result point will be projected along the segment to ensure 
                        // the result intersection point will not warp the segment
                        float2 point = segBeg.position + segBeg.direction * project;
                        dist = math.distancesq(vert.position, point);
                        if (dist < result.distanceSqr)
                        {
                            result.point = point;
                            result.distanceSqr = dist;
                            result.segment = segBeg.index;
                        }
                    }
                }
            }
        }
    }
}
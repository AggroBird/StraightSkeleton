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

        private sealed class Polygon
        {
            public Polygon(int index)
            {
                this.index = index;
                lhs = null;
                rhs = null;
                vertexCount = 0;
            }

            public readonly int index;

            // Bottom left starting vertex
            public PolygonVertex lhs;
            // Bottom right starting vertex
            public PolygonVertex rhs;
            // Total amount of vertices within this polygon
            public int vertexCount;


            public override string ToString()
            {
                return $"polygon #{index} ({vertexCount} vertices)";
            }
        }

        private sealed class PolygonVertex
        {
            public PolygonVertex(float2 position, int passIndex, Polygon polygon, PolygonSide side = PolygonSide.Unknown)
            {
                this.position = position;
                this.passIndex = passIndex;

                this.polygon = polygon;
                index = polygon.vertexCount++ | ((int)side << 30);
            }

            public Polygon polygon;
            private int index;

            public int Index
            {
                get => index & 0x3FFFFFFF;
                set => index = value | (index & ~0x3FFFFFFF);
            }
            public PolygonSide Side
            {
                get => (PolygonSide)((uint)index >> 30);
                set => index = Index | ((int)value << 30);
            }

            // The index of the pass this vertex was inserted
            public int passIndex;

            // 2D position of the vertex
            public float2 position;

            // Neighbours in the polygon linked list
            public PolygonVertex prevPolyVert;
            public PolygonVertex nextPolyVert;

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void Link(PolygonVertex prev, PolygonVertex next)
            {
                LinkPrev(prev);
                LinkNext(next);
            }
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void LinkPrev(PolygonVertex prev)
            {
                prevPolyVert = prev;
                prev.nextPolyVert = this;
            }
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void LinkNext(PolygonVertex next)
            {
                nextPolyVert = next;
                next.prevPolyVert = this;
            }


            public override string ToString()
            {
                switch (Side)
                {
                    case PolygonSide.Left:
                        return $"polygon #{polygon.index} vertex #{Index} (left side)";
                    case PolygonSide.Right:
                        return $"polygon #{polygon.index} vertex #{Index} (right side)";
                    default:
                        return $"polygon #{polygon.index} vertex #{Index}";
                }
            }
        }

        private sealed class ChainVertex
        {
            public ChainVertex(float2 position, int passIndex, int vertexIndex)
            {
                this.position = position;
                this.passIndex = passIndex;
                this.vertexIndex = vertexIndex;
            }

            public int vertexIndex;

            // The index of the pass this vertex was inserted
            public int passIndex;

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
            public ChainVertex prevChainVert;
            public ChainVertex nextChainVert;

            // Neighbouring polygons on this axis
            public PolygonVertex lhsPolyVert;
            public PolygonVertex rhsPolyVert;

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void Link(ChainVertex prev, ChainVertex next)
            {
                LinkPrev(prev);
                LinkNext(next);
            }
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void LinkPrev(ChainVertex prev)
            {
                prevChainVert = prev;
                prev.nextChainVert = this;
            }
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public void LinkNext(ChainVertex next)
            {
                nextChainVert = next;
                next.prevChainVert = this;
            }


            // Recalculate direction and length
            public void RecalculateSegment()
            {
                direction = nextChainVert.position - position;
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
            public void RecalculateBisector()
            {
                float2 a = direction;
                float2 b = -prevChainVert.direction;

                float signedAngle = math.atan2(b.x * a.y - b.y * a.x, math.dot(b, a));
                float PI2 = math.PI * 2;
                angle = (signedAngle + PI2) % PI2;
                float angleHalf = angle * 0.5f;
                bisector = Rotate(a, angleHalf);
                velocity = math.sin(angleHalf);
                IsReflex = angle > math.PI;
            }

            // Recalculate this node and its neighbours
            public void RecalculateNeighbours()
            {
                RecalculateSegment();
                prevChainVert.RecalculateSegment();
                nextChainVert.RecalculateSegment();
                RecalculateBisector();
                prevChainVert.RecalculateBisector();
                nextChainVert.RecalculateBisector();
            }


            public override string ToString()
            {
                return $"chain vertex #{vertexIndex}";
            }
        }

        // Chain of vertices.
        // Vertices are not guaranteed to be in order, so they are stored as a linked list.
        private class Chain : IEnumerable<ChainVertex>
        {
            public Chain()
            {

            }
            public Chain(IReadOnlyList<float2> points)
            {
                int pointCount = points.Count;
                for (int i = 0; i < pointCount; i++)
                {
                    vertices.Add(new ChainVertex(points[i], 0, i));
                }

                // Initialize the chain linked list
                for (int i = 0; i < pointCount; i++)
                {
                    int n = (i + 1) % pointCount;
                    int p = (i - 1 + pointCount) % pointCount;
                    vertices[i].nextChainVert = vertices[n];
                    vertices[i].prevChainVert = vertices[p];
                    vertices[i].RecalculateSegment();
                }
                // Initialize the bisectors
                for (int i = 0; i < pointCount; i++)
                {
                    vertices[i].RecalculateBisector();
                }
            }


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


            // Calculates the shortest shrink distance at which point a collision will occur
            public bool CalculateShortestShrinkDistance(out float shortest)
            {
                shortest = float.PositiveInfinity;
                for (int i = 0; i < vertices.Count; i++)
                {
                    ChainVertex v0 = vertices[i], v1 = v0.nextChainVert;

                    // Project the two bisectors onto each other and find the shortest collapse distance (edge event)
                    float2 n0 = new float2(-v0.bisector.y, v0.bisector.x);
                    float2 n1 = new float2(v1.bisector.y, -v1.bisector.x);

                    float entry0 = math.dot(v1.position - v0.position, n1) / math.dot(v0.bisector, n1);
                    float entry1 = math.dot(v0.position - v1.position, n0) / math.dot(v1.bisector, n0);

                    float entry = math.min(entry0 * v0.velocity, entry1 * v1.velocity);
                    if (!float.IsInfinity(entry) && entry > 0)
                    {
                        if (entry < shortest)
                        {
                            shortest = entry;
                        }
                    }

                    if (v0.IsReflex)
                    {
                        // Find the earliest segment intersect with a reflex vertex (split event)
                        v1 = v0;
                        do
                        {
                            ChainVertex v2 = v1.nextChainVert;
                            if (v1 != v0 && v2 != v0)
                            {
                                float2 perp = new float2(v1.direction.y, -v1.direction.x);
                                if (TracePlane(v1.position, perp, v0.position, v0.bisector, out entry))
                                {
                                    // Calculate time of collision (magic)
                                    float sin = math.sin(AngleBetween(v0.bisector, v1.direction));
                                    entry /= (1 / sin + 1 / v0.velocity);
                                    if (!float.IsInfinity(entry) && entry > 0)
                                    {
                                        // Grow the segment along the two bisectors of its vertices and project
                                        // the collision point to make sure this vertex will actually collide with
                                        // this segment in the future.
                                        float2 s0 = v1.position + v1.bisector * (entry / v1.velocity);
                                        float2 s1 = v2.position + v2.bisector * (entry / v2.velocity);
                                        float length = math.dot(s1 - s0, v1.direction);
                                        if (length > 0)
                                        {
                                            float2 point = v0.position + v0.bisector * (entry / v0.velocity);
                                            float project = math.dot(point - s0, v1.direction);
                                            if (project >= -MediumPrecisionEpsilon && project <= length + MediumPrecisionEpsilon)
                                            {
                                                if (entry < shortest)
                                                {
                                                    shortest = entry;
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            v1 = v2;
                        }
                        while (!v1.Equals(v0));
                    }
                }
                return shortest != float.PositiveInfinity;
            }

            // Shrink the polygon
            public void ApplyShrinkDistance(float distance)
            {
                for (int i = 0; i < vertices.Count; i++)
                {
                    ChainVertex vert = vertices[i];
                    vert.position += vert.bisector * (distance / vert.velocity);
                }

                // Recalculate the segments after shrinking
                for (int i = 0; i < vertices.Count; i++)
                {
                    vertices[i].RecalculateSegment();
                }
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
        }

        private int uniqueVertexCount = 0;
        private readonly List<Chain> activeChains = new List<Chain>();
        private readonly List<Chain> newChains = new List<Chain>();
        private readonly List<Polygon> polygons = new List<Polygon>();

        private readonly List<ChainVertex> incidentVertices = new List<ChainVertex>();
        private readonly List<ChainVertex> unresolvedChains = new List<ChainVertex>();
        private readonly List<ChainVertex> resolvedChains = new List<ChainVertex>();


        // Allow reusage of straight skeleton buffer data
        public void Generate(IReadOnlyList<float2> points, StraightSkeleton output)
        {
            output.Clear();

            if (points == null || points.Count < 3)
            {
                throw new StraightSkeletonException("Invalid input polygon provided");
            }

            Chain firstChain = new Chain(points);
            int pointCount = uniqueVertexCount = points.Count;

            polygons.Clear();
            for (int i = 0; i < pointCount; i++)
            {
                ChainVertex next = firstChain[i].nextChainVert;
                ChainVertex prev = next.prevChainVert;
                Polygon polygon = new Polygon(i);
                polygon.vertexCount = 0;
                polygon.lhs = new PolygonVertex(next.position, 0, polygon, PolygonSide.Left);
                polygon.rhs = new PolygonVertex(prev.position, 0, polygon, PolygonSide.Right);
                next.lhsPolyVert = polygon.lhs;
                prev.rhsPolyVert = polygon.rhs;
                polygon.lhs.prevPolyVert = polygon.rhs;
                polygon.rhs.nextPolyVert = polygon.lhs;
                polygons.Add(polygon);
            }

#if WITH_DEBUG
            foreach (var vert in firstChain)
            {
                output.AddDebugLine(vert.position, vert.nextChainVert.position, Color.white);
            }
#endif

            activeChains.Clear();
            activeChains.Add(firstChain);

            int passIndex = 0;
            while (activeChains.Count > 0)
            {
                if (passIndex++ > pointCount)
                {
                    throw new StraightSkeletonException("Max polygon shrink iteration count reached");
                }

                newChains.Clear();
                for (int i = 0; i < activeChains.Count; i++)
                {
                    Chain activeChain = activeChains[i];

                    if (!activeChain.CalculateShortestShrinkDistance(out float distance))
                    {
                        throw new StraightSkeletonException("Failed to find distance shortest shrink distance");
                    }

                    activeChain.ApplyShrinkDistance(distance);

                    ProcessIntersectionEvents(activeChain, passIndex);
                }

                activeChains.Clear();
                activeChains.AddRange(newChains);

#if WITH_DEBUG
                foreach (var chain in activeChains)
                {
                    foreach (var vert in chain)
                    {
                        output.AddDebugLine(vert.position, vert.nextChainVert.position, Color.grey);
                    }
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

                Polygon polygon = polygons[i];
                PolygonVertex first = polygon.lhs;
                PolygonVertex current = first;
                int iter = 0;
                bool hasInfiniteLoop = false;

                // Try and find the initial vertex (in case of unconnected loops)
                do
                {
                    if (iter++ > polygon.vertexCount)
                    {
                        Debug.LogError($"polygon #{polygon.index} contains an infinite loop");
                        hasInfiniteLoop = true;
                        isValid = false;
                        break;
                    }
                    if (current.prevPolyVert == null)
                    {
                        break;
                    }
                    current = current.prevPolyVert;
                }
                while (!current.Equals(first));
                first = current;

                // Fetch all vertices
                iter = 0;
                do
                {
                    if (iter++ > polygon.vertexCount)
                    {
                        if (!hasInfiniteLoop)
                        {
                            Debug.LogError($"polygon #{polygon.index} contains an infinite loop");
                            hasInfiniteLoop = true;
                        }
                        isValid = false;
                        break;
                    }
                    verts.Add(current.position);
                    if (current.nextPolyVert == null)
                    {
                        Debug.LogError($"polygon #{polygon.index} contains null connections");
                        isValid = false;
                        break;
                    }
                    current = current.nextPolyVert;
                }
                while (!current.Equals(first));

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
                PolygonVertex first = polygons[i].lhs;
                PolygonVertex vertex = first;
                do
                {
                    buffer.Add(vertex.position);
                    vertex = vertex.nextPolyVert;
                }
                while (!vertex.Equals(first));
            }
        }
        public StraightSkeleton Generate(IReadOnlyList<float2> points)
        {
            StraightSkeleton output = new StraightSkeleton();
            Generate(points, output);
            return output;
        }

        // Check for any intersections (results in one or more new chains)
        private void ProcessIntersectionEvents(Chain chain, int passIndex)
        {
            ChainVertex vertex;

            // This phase scans over all the vertices and checks for overlapping (incident) vertices.
            // In the case of a split event, it will emit a new vertex parallel to the segment so that
            // it can be picked up later by the split stage.
            ChainVertex first = chain.First;
            vertex = first;
            do
            {
                // No need to resolve intersections with inserted split polygons
                if (vertex.passIndex != passIndex)
                {
                    // Find closest intersect
                    IntersectResult bestResult = IntersectResult.None;
                    ChainVertex segment = vertex;
                    do
                    {
                        if (segment != vertex && segment.nextChainVert != vertex)
                        {
                            AppendIntersect(vertex, segment, LowPrecisionEpsilon, ref bestResult);
                        }

                        segment = segment.nextChainVert;
                    }
                    while (!segment.Equals(vertex));

                    if (bestResult)
                    {
                        // Snap vertex to intersect location
                        vertex.position = bestResult.point;
                        vertex.RecalculateNeighbours();

                        if (!bestResult.IsIncident)
                        {
                            // Insert new vertex on segment
                            ChainVertex prev = bestResult.segment, next = prev.nextChainVert;
                            ChainVertex insert = new ChainVertex(bestResult.point, passIndex, uniqueVertexCount++);
                            // Create an empty polygon vertex at the split position
                            PolygonVertex polygonVertex = new PolygonVertex(bestResult.point, passIndex, prev.rhsPolyVert.polygon);
                            insert.lhsPolyVert = insert.rhsPolyVert = polygonVertex;
                            insert.Link(prev, next);
                            insert.RecalculateNeighbours();
                        }
                    }
                }

                vertex = vertex.nextChainVert;
            }
            while (!vertex.Equals(first));

            // This stage takes all incident vertices and splits them up into seperate chains.
            // It also emits polygon vertices for the skeleton in the case of an edge or split event.
            resolvedChains.Clear();
            unresolvedChains.Clear();
            unresolvedChains.Add(chain.First);
        TraceNext:
            while (unresolvedChains.Count > 0)
            {
                int index = unresolvedChains.Count - 1;
                ChainVertex unresolved = unresolvedChains[index];
                unresolvedChains.RemoveAt(index);

                vertex = unresolved;
                do
                {
                    // No need to resolve vertices that we have already checked for incident in this pass
                    if (vertex.passIndex != passIndex)
                    {
                        vertex.passIndex = passIndex;

                        incidentVertices.Clear();

                        // Collect all vertices that are incident to the current vertex
                        ChainVertex other = vertex;
                        do
                        {
                            if (vertex != other && IncidentVertices(vertex, other, LowPrecisionEpsilon * 2))
                            {
                                other.passIndex = passIndex;

                                incidentVertices.Add(other);
                            }

                            other = other.nextChainVert;
                        }
                        while (!other.Equals(vertex));

                        if (incidentVertices.Count > 0)
                        {
                            // Resolve the current vertex as well
                            incidentVertices.Add(vertex);

                            // Split into subshapes
                            int incidentCount = incidentVertices.Count;
                            for (int i = 0; i < incidentVertices.Count; i++)
                            {
                                ChainVertex prevIncident = incidentVertices[i];
                                ChainVertex nextIncident = incidentVertices[(i + 1) % incidentCount];

                                // If there are vertices between the previous and next incident vertices,
                                // we have to split them off into a new chain
                                if (prevIncident.nextChainVert != nextIncident)
                                {
                                    ChainVertex insert = new ChainVertex(vertex.position, passIndex, uniqueVertexCount++);

                                    PolygonVertex nextPolygon = prevIncident.rhsPolyVert;
                                    if (nextPolygon.passIndex != passIndex)
                                    {
                                        // Insert a polygon vertex on the right side of the left polygon
                                        PolygonVertex lhs = new PolygonVertex(vertex.position, passIndex, nextPolygon.polygon);
                                        lhs.LinkNext(nextPolygon);
                                        insert.rhsPolyVert = lhs;
                                    }
                                    else
                                    {
                                        // If the polygon point inserted through a split event, we dont have to insert a new vertex
                                        insert.rhsPolyVert = nextPolygon;
                                    }

                                    PolygonVertex prevPolygon = nextIncident.lhsPolyVert;
                                    if (prevPolygon.passIndex != passIndex)
                                    {
                                        // Insert a polygon vertex on the left side of the right polygon
                                        PolygonVertex rhs = new PolygonVertex(vertex.position, passIndex, prevPolygon.polygon);
                                        rhs.LinkPrev(prevPolygon);
                                        insert.lhsPolyVert = rhs;
                                    }
                                    else
                                    {
                                        // If the polygon point inserted through a split event, we dont have to insert a new vertex
                                        insert.lhsPolyVert = prevPolygon;
                                    }

                                    insert.Link(nextIncident.prevChainVert, prevIncident.nextChainVert);
                                    insert.RecalculateNeighbours();
                                    unresolvedChains.Add(insert);
                                }
                                else
                                {
                                    // Close into a triangle
                                    PolygonVertex end = new PolygonVertex(vertex.position, passIndex, prevIncident.rhsPolyVert.polygon);
                                    end.Link(nextIncident.lhsPolyVert, prevIncident.rhsPolyVert);
                                }
                            }

                            goto TraceNext;
                        }
                    }

                    vertex = vertex.nextChainVert;
                }
                while (!vertex.Equals(unresolved));

                resolvedChains.Add(unresolved);
            }

            // Final stage: rebuild new chains from the split results
            for (int i = 0; i < resolvedChains.Count; i++)
            {
                Chain newChain = new Chain();
                first = resolvedChains[i];
                vertex = first;
                do
                {
                    newChain.Add(vertex);

                    vertex = vertex.nextChainVert;
                }
                while (!vertex.Equals(first));

                if (newChain.Count > 2)
                {
                    // New valid chain
                    newChains.Add(newChain);
                }
                else if (newChain.Count == 2)
                {
                    // If we have only two vertices we can collapse the segment and close the polygons
                    ChainVertex a = newChain[0], b = newChain[1];
                    if (a.lhsPolyVert.nextPolyVert == null && b.rhsPolyVert.prevPolyVert == null)
                    {
                        a.lhsPolyVert.LinkNext(b.rhsPolyVert);
                    }
                    if (a.rhsPolyVert.prevPolyVert == null && b.lhsPolyVert.nextPolyVert == null)
                    {
                        a.rhsPolyVert.LinkPrev(b.lhsPolyVert);
                    }
                }
            }

            chain.Clear();
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
        private static bool IncidentVertices(ChainVertex a, ChainVertex b, float accuracy)
        {
            return IncidentVertices(a.position, b.position, accuracy);
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
            public static readonly IntersectResult None = new IntersectResult { distanceSqr = float.PositiveInfinity };

            public float2 point;
            public float distanceSqr;
            public ChainVertex segment;

            public bool IsIncident => segment == null;

            public static implicit operator bool(IntersectResult intersectResult)
            {
                return intersectResult.distanceSqr != float.PositiveInfinity;
            }
        }

        // Calculate the intersection of a vertex on a segment, and if the distance is closer than the previous, update the result.
        // This function prefers incident results, to reduce the amount of tiny segment intersections.
        private static void AppendIntersect(ChainVertex vertex, ChainVertex segment, float accuracy, ref IntersectResult result)
        {
            float acc2 = accuracy * accuracy;
            ChainVertex v2 = segment.nextChainVert;
            float2 p0 = vertex.position;
            float2 p1 = segment.position;
            float2 p2 = v2.position;

            // See if we are anywhere near the segment origin
            float dist = math.distancesq(p0, p1);
            if (dist <= acc2 && dist < result.distanceSqr)
            {
                result.point = p1;
                result.distanceSqr = dist;
                result.segment = null;
            }

            // See if we are anywhere near the segment end
            dist = math.distancesq(p0, p2);
            if (dist <= acc2 && dist < result.distanceSqr)
            {
                result.point = p2;
                result.distanceSqr = dist;
                result.segment = null;
            }

            // Project the point onto the segment (if we are currently not incident)
            if (segment.length > 0 && (!result || !result.IsIncident))
            {
                // Project on the segment along the length
                float2 relative = vertex.position - segment.position;
                float project = math.dot(relative, segment.direction);
                if (project >= 0 && project <= segment.length)
                {
                    // Project on the segment along the width (with an error margin)
                    float2 perp = new float2(segment.direction.y, -segment.direction.x);
                    float dot = math.dot(perp, relative);
                    if (math.abs(dot) <= accuracy)
                    {
                        // The result point will be projected along the segment to ensure 
                        // the result intersection point will not warp the segment
                        float2 point = segment.position + segment.direction * project;
                        dist = math.distancesq(vertex.position, point);
                        if (dist < result.distanceSqr)
                        {
                            result.point = point;
                            result.distanceSqr = dist;
                            result.segment = segment;
                        }
                    }
                }
            }
        }
    }
}
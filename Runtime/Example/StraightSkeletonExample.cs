using System;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

namespace AggroBird.StraightSkeleton
{
    [ExecuteInEditMode, RequireComponent(typeof(MeshRenderer), typeof(MeshFilter))]
    public class StraightSkeletonExample : MonoBehaviour
    {
        public bool snapToGrid = false;
        public Vector2[] polygon = new Vector2[0];
        [Space]
        [Min(0)] public float wallHeight = 2;
        public Color wallColor = Color.white;
        [Range(1, 89)] public float roofPitch = 45;
        public Color roofColor = Color.red;
        [Space]
        public bool generateOnValidate = false;

        public float maxHeight = 999;


        private readonly StraightSkeletonGenerator generator = new StraightSkeletonGenerator();
        private readonly StraightSkeleton straightSkeleton = new StraightSkeleton();

        private readonly List<Vector3> vertices = new List<Vector3>();
        private readonly List<Vector2> texcoords = new List<Vector2>();
        private readonly List<Color> colors = new List<Color>();
        private readonly List<int> triangles = new List<int>();

        private struct TriangulationNode
        {
            public TriangulationNode(int idx, float2 coord, int prev, int next)
            {
                this.idx = idx;
                this.coord = coord;
                this.prev = prev;
                this.next = next;
            }

            public int idx;
            public float2 coord;
            public int prev;
            public int next;
        }
        private TriangulationNode[] triangulation = new TriangulationNode[32];

        private bool updateMesh = false;


        public void GenerateBuilding()
        {
            if (polygon != null && polygon.Length > 2)
            {
                float2[] input = new float2[polygon.Length];
                for (int i = 0; i < polygon.Length; i++)
                {
                    input[i] = polygon[i];
                }

                // Generate straight skeleton
                try
                {
                    generator.Generate(input, straightSkeleton, maxHeight);
                }
                catch (Exception ex)
                {
                    Debug.LogException(ex);
                    return;
                }

                GenerateMesh();
            }
        }

        private void GenerateMesh()
        {
            MeshFilter meshFilter = GetComponent<MeshFilter>();
            if (!meshFilter)
            {
                Debug.LogError("Failed to find mesh filter component");
                return;
            }
            Mesh mesh = meshFilter.sharedMesh ? meshFilter.sharedMesh : new Mesh();
            mesh.hideFlags |= HideFlags.DontSave | HideFlags.NotEditable;

            vertices.Clear();
            texcoords.Clear();
            colors.Clear();
            triangles.Clear();

            // Generate walls
            if (wallHeight > 0)
            {
                for (int i = 0, k = 0; i < polygon.Length; i++, k += 4)
                {
                    int j = (i + 1) % polygon.Length;
                    Vector2 p0 = polygon[i];
                    Vector2 p1 = polygon[j];
                    float size = math.length(p1 - p0);

                    vertices.Add(new Vector3(p0.x, 0, p0.y));
                    vertices.Add(new Vector3(p0.x, wallHeight, p0.y));
                    vertices.Add(new Vector3(p1.x, wallHeight, p1.y));
                    vertices.Add(new Vector3(p1.x, 0, p1.y));
                    texcoords.Add(new Vector2(size, 0));
                    texcoords.Add(new Vector2(size, wallHeight));
                    texcoords.Add(new Vector2(0, wallHeight));
                    texcoords.Add(new Vector2(0, 0));
                    colors.Add(wallColor);
                    colors.Add(wallColor);
                    colors.Add(wallColor);
                    colors.Add(wallColor);
                    triangles.Add(k + 0);
                    triangles.Add(k + 2);
                    triangles.Add(k + 1);
                    triangles.Add(k + 3);
                    triangles.Add(k + 2);
                    triangles.Add(k + 0);
                }
            }

            // Generate roof
            float slope = math.tan(math.radians(roofPitch));
            float v = 1.0f / math.cos(math.radians(roofPitch));
            for (int idx = 0; idx < straightSkeleton.Count; idx++)
            {
                var subPolygon = straightSkeleton[idx];
                float texcoordScale = idx < straightSkeleton.CenterPolygonCount ? 1 : v;
                int vertCount = subPolygon.Count;
                if (vertCount < 3) continue;

                if (vertCount > triangulation.Length)
                {
                    triangulation = new TriangulationNode[vertCount];
                }
                float2 origin = subPolygon[0].xy;
                float2 dir = math.normalize(origin - subPolygon[vertCount - 1].xy);
                float2 perp = new float2(dir.y, -dir.x);
                int offset = vertices.Count;
                for (int i = 0; i < vertCount; i++)
                {
                    int next = (i + 1) % vertCount;
                    int prev = (i + vertCount - 1) % vertCount;
                    float3 vert = subPolygon[i];
                    float2 pos = vert.xy;
                    float2 relative = pos - origin;
                    triangulation[i] = new TriangulationNode(offset + i, pos, prev, next);
                    vertices.Add(new Vector3(pos.x, vert.z * slope + wallHeight, pos.y));
                    texcoords.Add(new Vector2(math.dot(relative, dir), math.dot(relative, perp) * texcoordScale));
                    colors.Add(roofColor);
                }

                // Quick and dirty ear clipping algorithm
                int remainingCount = vertCount;
                int outer = 0;
            ProcessNext:
                if (remainingCount > 2)
                {
                    for (int i = 0; i < remainingCount; i++)
                    {
                        TriangulationNode node = triangulation[outer];
                        TriangulationNode next = triangulation[node.next];
                        TriangulationNode prev = triangulation[node.prev];
                        float2 c0 = node.coord;
                        float2 c1 = next.coord;
                        float2 c2 = prev.coord;
                        float2 d0 = c1 - c0;
                        float2 d2 = c0 - c2;
                        float2 p2 = new float2(d2.y, -d2.x);
                        if (math.dot(p2, d0) > 0)
                        {
                            if (remainingCount > 3)
                            {
                                float2 d1 = c2 - c1;
                                float2 p0 = new float2(d0.y, -d0.x);
                                float2 p1 = new float2(d1.y, -d1.x);
                                int otherCount = remainingCount - 3;
                                int inner = next.next;
                                for (int j = 0; j < otherCount; j++)
                                {
                                    float2 p = triangulation[inner].coord;
                                    if (math.dot(p - c0, p0) > 0)
                                    {
                                        if (math.dot(p - c1, p1) > 0)
                                        {
                                            if (math.dot(p - c2, p2) > 0)
                                            {
                                                goto Skip;
                                            }
                                        }
                                    }
                                    inner = triangulation[inner].next;
                                }
                            }

                            triangles.Add(node.idx);
                            triangles.Add(next.idx);
                            triangles.Add(prev.idx);

                            triangulation[node.next].prev = node.prev;
                            triangulation[node.prev].next = node.next;
                            outer = node.next;

                            remainingCount--;
                            goto ProcessNext;
                        }

                    Skip:
                        outer = node.next;
                    }
                }
            }

            // Assign mesh data
            mesh.Clear();
            mesh.SetVertices(vertices);
            mesh.SetUVs(0, texcoords);
            mesh.SetColors(colors);
            mesh.SetTriangles(triangles, 0);
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();
            meshFilter.sharedMesh = mesh;
        }


        private void Start()
        {
            GenerateBuilding();
        }

        private void Update()
        {
            if (updateMesh)
            {
                updateMesh = false;

                GenerateBuilding();
            }
        }

        private void OnValidate()
        {
            if (generateOnValidate)
            {
                // Not allowed to assign mesh in OnValidate
                updateMesh = true;
            }
        }

        private void OnDrawGizmosSelected()
        {
            float2 scale = new float2(transform.lossyScale.x, transform.lossyScale.z);
            foreach (var line in straightSkeleton.DebugOutput)
            {
                Gizmos.color = line.color;
                Vector2 p0 = line.from * scale, p1 = line.to * scale;
                Gizmos.DrawLine(new Vector3(p0.x, 0, p0.y), new Vector3(p1.x, 0, p1.y));
            }
        }
    }
}
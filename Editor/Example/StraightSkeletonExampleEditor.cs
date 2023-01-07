using UnityEditor;
using UnityEngine;

namespace AggroBird.StraightSkeleton.Editor
{
    [CustomEditor(typeof(StraightSkeletonExample))]
    public class StraightSkeletonExampleEditor : UnityEditor.Editor
    {
        public override void OnInspectorGUI()
        {
            StraightSkeletonExample example = target as StraightSkeletonExample;
            serializedObject.Update();
            SerializedProperty polygon = serializedObject.FindProperty("polygon");
            if (GUILayout.Button("Generate"))
            {
                example.GenerateBuilding();
                SceneView.RepaintAll();
            }
            if (GUILayout.Button("Add Point"))
            {
                polygon.arraySize++;
                if (polygon.arraySize >= 3)
                {
                    Vector2 v0 = polygon.GetArrayElementAtIndex(polygon.arraySize - 3).vector2Value;
                    Vector2 v1 = polygon.GetArrayElementAtIndex(polygon.arraySize - 2).vector2Value;
                    polygon.GetArrayElementAtIndex(polygon.arraySize - 1).vector2Value = v1 + (v1 - v0).normalized;
                }
            }
            EditorGUILayout.BeginHorizontal();
            if (GUILayout.Button("Shift backward"))
            {
                Vector2[] copy = CopyArray();
                for (int i = 0; i < copy.Length; i++)
                {
                    polygon.GetArrayElementAtIndex((i + 1) % copy.Length).vector2Value = copy[i];
                }
            }
            if (GUILayout.Button("Shift forward"))
            {
                Vector2[] copy = CopyArray();
                for (int i = 0; i < copy.Length; i++)
                {
                    polygon.GetArrayElementAtIndex((i + copy.Length - 1) % copy.Length).vector2Value = copy[i];
                }
            }
            EditorGUILayout.EndHorizontal();
            serializedObject.ApplyModifiedProperties();
            base.OnInspectorGUI();
        }
        private Vector2[] CopyArray()
        {
            StraightSkeletonExample example = target as StraightSkeletonExample;
            if (example.polygon == null || example.polygon.Length == 0) return System.Array.Empty<Vector2>();
            Vector2[] result = new Vector2[example.polygon.Length];
            System.Array.Copy(example.polygon, result, example.polygon.Length);
            return result;
        }

        public void OnSceneGUI()
        {
            Camera camera = SceneView.currentDrawingSceneView.camera;
            StraightSkeletonExample example = target as StraightSkeletonExample;
            if (camera && example)
            {
                serializedObject.Update();
                SerializedProperty polygon = serializedObject.FindProperty("polygon");
                Handles.color = Color.white;
                Vector3 origin = camera.transform.position;
                Vector3 labelOffset = camera.transform.right - camera.transform.up;
                bool snapToGrid = example.snapToGrid;
                for (int i = 0; i < polygon.arraySize; i++)
                {
                    SerializedProperty value = polygon.GetArrayElementAtIndex(i);
                    Vector2 val2D = value.vector2Value;
                    Vector3 currentPosition = example.transform.TransformPoint(new Vector3(val2D.x, 0, val2D.y));
                    float scale = Vector3.Distance(origin, currentPosition) * 0.02f;
                    Vector3 newPosition = Handles.FreeMoveHandle(currentPosition, Quaternion.identity, scale, Vector3.zero, Handles.CircleHandleCap);
                    if (currentPosition != newPosition)
                    {
                        Vector3 direction = (newPosition - origin).normalized;
                        float entry = Vector3.Dot(example.transform.position - origin, example.transform.up) / Vector3.Dot(direction, example.transform.up);
                        if (entry > 0)
                        {
                            newPosition = example.transform.InverseTransformPoint(origin + direction * entry);
                            val2D = new Vector2(newPosition.x, newPosition.z);
                            if (snapToGrid)
                            {
                                val2D.x = Mathf.Round(val2D.x);
                                val2D.y = Mathf.Round(val2D.y);
                            }
                            value.vector2Value = val2D;
                        }
                    }
                    newPosition += labelOffset * scale;
                    Handles.Label(newPosition, i.ToString());
                }
                serializedObject.ApplyModifiedProperties();
            }
        }
    }
}
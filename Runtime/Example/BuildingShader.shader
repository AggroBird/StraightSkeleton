Shader "Custom/BuildingShader"
{
    Properties
    {
        [NoScaleOffset] _MainTex("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 200

        CGPROGRAM
        #pragma surface surf Standard fullforwardshadows vertex:vert
        #pragma target 3.0

        sampler2D _MainTex;

        struct Input
        {
            float2 uv_MainTex;
            float4 color : COLOR;
        };

        void vert(inout appdata_full v, out Input o)
        {
            UNITY_INITIALIZE_OUTPUT(Input, o);
        }

        void surf(Input input, inout SurfaceOutputStandard output)
        {
            output.Albedo = tex2D(_MainTex, input.uv_MainTex).rgb * input.color.rgb;
            output.Alpha = 1;
        }
        ENDCG
    }
    FallBack "Diffuse"
}

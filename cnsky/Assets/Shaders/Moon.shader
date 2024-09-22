Shader "cnlohrdomnomnom/Moon"
{
	//based on Shader "d4rkpl4y3r/BRDF PBS Macro"
	
	Properties
	{
		[Enum(Off, 0, Front, 1, Back, 2)] _Culling ("Culling Mode", Int) = 2
		_Cutoff("Cutout", Range(0,1)) = .5
		_MainTexDay("DayText", 2D) = "white" {}
		_MainTexNight("NightTex", 2D) = "white" {}
		[hdr] _Color("Albedo", Color) = (1,1,1,1)
		[Gamma] _Metallic("Metallic", Range(0, 1)) = 0
		_Smoothness("Smoothness", Range(0, 1)) = 0
		_EarthMidnightRotation("Midnight Rotation", float) = 0
		_ManagementTexture ("Management", 2D) = "white" {}
		[UIToggle] _RayTrace ("Ray Trace", float) = 1.0
	}
	SubShader
	{
		Tags
		{
			"RenderType"="Opaque"
			"Queue"="Geometry"
		}

		Cull [_Culling]

		CGINCLUDE
			#include "UnityCG.cginc"
			#include "Lighting.cginc"
			#include "AutoLight.cginc"
			#include "UnityPBSLighting.cginc"
			#include "Packages/com.llealloo.audiolink/Runtime/Shaders/AudioLink.cginc"

			uniform float4 _Color;
			uniform float _Metallic;
			uniform float _Smoothness;
			float _RayTrace;
			uniform sampler2D _MainTexDay, _MainTexNight;
			//uniform sampler2D _MainTexDay, _MainTexNight;
			uniform float4 _MainTexDay_ST;
			uniform float _Cutoff, _EarthMidnightRotation;
			float4 _ManagementTexture_TexelSize;
			Texture2D< float4 > _ManagementTexture;

			struct v2f
			{
				float3 rayOrigin : RAYORIGIN;
				float3 rayDir : RAYDIR;

				#ifndef UNITY_PASS_SHADOWCASTER
				float4 pos : SV_POSITION;
				float3 normal : NORMAL;
				float3 wPos : WPOSC;
				float3 vPos : VPOSC;
				
				
				SHADOW_COORDS(3)
				#else
				V2F_SHADOW_CASTER;
				#endif
				float2 uv : TEXCOORD1;
			};


			// def sun_position(mjd_utc, dt_tt=69.184):
			//     """ECI position (meters) of the Sun at the given MJD (UTC)."""
			//     # Reference: Meeus, Jean, Astronomical Algorithms (2nd Ed.). Richmond:
			//     # Willmann-Bell, Inc., 2009.
			//     T = (mjd_utc - 51544.5 + dt_tt / 86400) / 36525                  # (25.1)
			//     L0_deg = 280.46646 + T * (36000.76983 + T * 0.0003032)           # (25.2)
			//     M_deg = 357.52911 + T * (35999.05029 - T * 0.0001537)            # (25.3)
			//     M = radians(M_deg)
			//     e = 0.016708634 - T * (0.000042037 + T * 0.0000001267)           # (25.4)
			//     C_deg = ((1.914602 - T * (0.004817 + T * 0.000014)) * sin(M)
			//            + (0.019993 - T * 0.000101) * sin(2 * M)
			//             + 0.000289 * sin(3 * M))
			//     true_lon_deg = L0_deg + C_deg
			//     nu = radians(M_deg + C_deg)
			//     R = 149597870000 * 1.000001018 * (1 - e * e) / (1 + e * cos(nu)) # (25.5)
			//     Omega = radians(125.04 - 1934.136 * T)
			//     lon = radians(true_lon_deg - 0.00569 - 0.00478 * sin(Omega))
			//     epsilon = radians(23 + (26 + (21.448                             # (22.2)
			//         - T * (46.8150 + T * (0.00059 - T * 0.001813))) / 60) / 60
			//         + 0.00256 * cos(Omega))                                      # (25.8)
			//     return (R * cos(lon),                                            # (26.1)
			//             R * sin(lon) * cos(epsilon),
			//             R * sin(lon) * sin(epsilon))

			float radians(float degrees) {
				return degrees * (2.*3.1415926535 / 360.);
			}

			float3 moon_position(float days_since_J2000) {

			    // ECI position (meters) of the Moon at the given MJD (UTC).
			    // ported from public domain code:
			    // https://possiblywrong.wordpress.com/2014/01/04/computing-positions-and-eclipses-of-the-sun-and-moon/
			    // Reference: Montenbruck and Gill, Satellite Orbits. Berlin: Springer,
			    // 2005, Chapter 3.3.2.
			    // float days_since_J2000 = mjd_utc - 51544.5 + dt_tt / 86400;
			    float T = days_since_J2000 / 36525.; // centuries since J2000

			    // equation (3.47)
			    float L0 = radians(218.31617 + 481267.88088 * T);
			    float l = radians(134.96292 + 477198.86753 * T);
			    float lp = radians(357.52543 + 35999.04944 * T);
			    float F = radians(93.27283 + 483202.01873 * T);
			    float D = radians(297.85027 + 445267.11135 * T);

			    // equation (3.48)
			    float dL = radians((22640 * sin(l) + 769 * sin(2 * l)
			        - 4586 * sin(l - 2 * D) + 2370 * sin(2 * D)
			        - 668 * sin(lp) - 412 * sin(2 * F)
			        - 212 * sin(2 * l - 2 * D) - 206 * sin(l + lp - 2 * D)
			        + 192 * sin(l + 2 * D) - 165 * sin(lp - 2 * D)
			        + 148 * sin(l - lp) - 125 * sin(D)
			        - 110 * sin(l + lp) - 55 * sin(2 * F - 2 * D)) / 3600);
			    float lon = L0 + dL;

			    // equation (3.49)
			    float beta = radians((18520 * sin(F + dL + radians(
			            (412 * sin(2 * F) + 541 * sin(lp)) / 3600))
			        - 526 * sin(F - 2 * D) + 44 * sin(l + F - 2 * D)
			        - 31 * sin(-l + F - 2 * D) - 25 * sin(-2 * l + F)
			        - 23 * sin(lp + F - 2 * D) + 21 * sin(-l + F)
			        + 11 * sin(-lp + F - 2 * D)) / 3600);

			    // equation (3.50)
			    float r = 1000 * (385000 - 20905 * cos(l) - 3699 * cos(2 * D - l)
			         - 2956 * cos(2 * D) - 570 * cos(2 * l) + 246 * cos(2 * l - 2 * D)
			         - 205 * cos(lp - 2 * D) - 171 * cos(l + 2 * D)
			         - 152 * cos(l + lp - 2 * D));

			    // equation (3.51)
			    float x = r * cos(lon) * cos(beta);
			    float y = r * sin(lon) * cos(beta);
			    float z = r * sin(beta);

			    // equation (3.45)
			    float epsilon = radians(23.43929111);
			    float s = sin(-epsilon);
			    float c = cos(-epsilon);
			    return float3(
			    	x,
			    	 y * c + z * s,
			    	-y * s + z * c
		    	);
			}

			v2f vert(appdata_base v)
			{
				v2f o;

				float UTCDAY = AudioLinkDecodeDataAsUInt( ALPASS_GENERALVU_UNIX_DAYS );
				float UTCDAYf = AudioLinkDecodeDataAsSeconds( ALPASS_GENERALVU_UNIX_SECONDS )/86400.0;
				float J2000_in_unix_days = (946684800 / 86400.0);
				float days_since_J2000 = (UTCDAY - J2000_in_unix_days) + UTCDAYf;

				// Move the moon to a new position
				float4x4 ObjectToWorld = unity_ObjectToWorld;
				float3 objectOrigin = mul(ObjectToWorld, float4(0.0,0.0,0.0,1.0) );
				float distance_to_earth = length(objectOrigin);
				float3 moon_eci_meters = moon_position(days_since_J2000);
				objectOrigin = distance_to_earth * normalize(moon_eci_meters); //* float3(sin(_Time.y), 0, cos(_Time.y));
				ObjectToWorld[0][3] = objectOrigin.x;
				ObjectToWorld[1][3] = objectOrigin.y;
				ObjectToWorld[2][3] = objectOrigin.z;

				// Thanks ben dot com.
				// I saw these ortho shadow substitutions in a few places, but bgolus explains them
				// https://bgolus.medium.com/rendering-a-sphere-on-a-quad-13c92025570c
				float howOrtho = UNITY_MATRIX_P._m33; // instead of unity_OrthoParams.w
				float3 worldSpaceCameraPos = UNITY_MATRIX_I_V._m03_m13_m23; // instead of _WorldSpaceCameraPos
				float3 worldPos = mul(ObjectToWorld, v.vertex);
				float3 cameraToVertex = worldPos - worldSpaceCameraPos;
				float3 orthoFwd = -UNITY_MATRIX_I_V._m02_m12_m22; // often seen: -UNITY_MATRIX_V[2].xyz;
				float3 orthoRayDir = orthoFwd * dot(cameraToVertex, orthoFwd);
				// start from the camera plane (can also just start from o.vertex if your scene is contained within the geometry)
				float3 orthoCameraPos = worldPos - orthoRayDir;

				o.rayOrigin = lerp(worldSpaceCameraPos, orthoCameraPos, howOrtho );
				o.rayDir = normalize( lerp( cameraToVertex, orthoRayDir, howOrtho ) );
				
				// Switch to view
				o.rayOrigin = mul( UNITY_MATRIX_V, float4( o.rayOrigin, 1.0 ) );
				o.rayDir = mul( UNITY_MATRIX_V, float4( o.rayDir, 0.0 ) );
				
				// Switch back to object.
				o.rayOrigin = mul( float4( o.rayOrigin, 1.0 ), UNITY_MATRIX_IT_MV );
				o.rayDir = mul( float4( o.rayDir, 0.0 ), UNITY_MATRIX_IT_MV ); //?!?!?! Why broke?
				
				
				#ifdef UNITY_PASS_SHADOWCASTER
				TRANSFER_SHADOW_CASTER_NOPOS(o, o.pos);
				#else
				o.wPos = mul(ObjectToWorld, v.vertex);
				o.vPos = v.vertex;
				o.pos = UnityWorldToClipPos(o.wPos);
				o.normal = UnityObjectToWorldNormal(v.normal);
				TRANSFER_SHADOW(o);
				#endif
				o.uv = TRANSFORM_TEX(v.texcoord.xy, _MainTexDay);
				return o;
			}

			#ifndef UNITY_PASS_SHADOWCASTER
			float4 frag(v2f i) : SV_TARGET
			{
				float3 normal = normalize(i.normal);
				float2 uv = i.uv;
				float2 uvbase = uv;
				float3 wPos = i.wPos.xyz;
				float3 vPos = i.vPos.xyz;
				
				float4 InfoBlock0 = _ManagementTexture.Load( int3( 0, _ManagementTexture_TexelSize.w - 1, 0 ) );
				float fFTime = InfoBlock0.z;
				
				float lambda;
				float phi;
				
				if( _RayTrace )
				{
					float3 rayOrigin = i.rayOrigin;
					float3 rayDir = normalize( vPos - rayOrigin );//(i.rayDir);

					float4 sphere = float4( 0, 0, 0, 0.5 );
					float bestHitDistance = 1e20;
					float3 bestHitNormal = 0.0;
					float3 bestHitPosition = 0.0;
					// Calculate distance along the ray where the sphere is intersected
					float3 d = rayOrigin - sphere.xyz;
					float p1 = -dot(rayDir, d);
					float p2sqr = p1 * p1 - dot(d, d) + sphere.w * sphere.w;
					if (p2sqr < 0)
						discard;
					float p2 = sqrt(p2sqr);
					float t = p1 - p2 > 0 ? p1 - p2 : p1 + p2;
					if (t > 0 && t < bestHitDistance)
					{
						bestHitPosition = rayOrigin + t * rayDir;
						bestHitNormal = normalize(bestHitPosition - sphere.xyz);
						bestHitDistance = t;
					}
					else
					{
						discard;
					}

	
					lambda = atan2( bestHitPosition.z, bestHitPosition.x );
					phi = atan2( length(bestHitPosition.xz), bestHitPosition.y );
					vPos = bestHitPosition;
					wPos = mul( UNITY_MATRIX_M, float4( bestHitPosition, 1.0 ) );
					normal = normalize( mul( UNITY_MATRIX_M, float4( bestHitNormal, 0.0 ) ) );
				}
				else
				{
					lambda = atan2( i.vPos.z, i.vPos.x );
					phi = atan2( length(i.vPos.xz), i.vPos.y );
				}
				
				lambda -= fFTime*3.1415926535*2 + _EarthMidnightRotation;

				uv.x = (lambda/2.0);
				uv.y = -phi;

				uv.x = frac( uv.x / 3.1415926535 + 1 );
				uv.y = frac( uv.y / 3.1415926535 + 1 );
				//uv = clamp( uv, 0.0, .99 );


				float3 sun_dir;// = normalize(_WorldSpaceLightPos0.xyz);
				{   // TODO: do this once. maybe in udon or in the management texture.
					// https://astronomy.stackexchange.com/a/37199
					float UTCDAY = AudioLinkDecodeDataAsUInt( ALPASS_GENERALVU_UNIX_DAYS );
					float UTCDAYf = AudioLinkDecodeDataAsSeconds( ALPASS_GENERALVU_UNIX_SECONDS )/86400.0;
					float J2000_in_unix_days = (946684800 / 86400.0);
					float d = (UTCDAY - J2000_in_unix_days) + UTCDAYf;  // days_since_J2000

					const float pi = 3.14159265359;
				    float L = 280.4606184 + ((36000.77005361 / 36525) * d); // mean longitude, in degrees
				    float g = 357.5277233 + ((35999.05034 / 36525) * d); // mean anomaly, in degrees
				    float p = L + (1.914666471 * sin(g * pi / 180)) + (0.918994643 * sin(2*g * pi / 180)); // ecliptic longitude lambda, in degrees
				    float q = 23.43929 - ((46.8093/3600) * (d / 36525)); // obliquity of ecliptic plane epsilon, in degrees

				    sun_dir = float3(
					    cos(p * pi / 180),
					    cos(q * pi / 180) * sin(p * pi / 180),
					    sin(q * pi / 180) * sin(p * pi / 180)
					);
				}
				float dayness = saturate( dot( normal,  sun_dir)  * 6.0);
				float4 texCol = lerp( tex2Dgrad(_MainTexNight, uv, ddx(uvbase), ddy(uvbase) ), tex2Dgrad(_MainTexDay, uv, ddx(uvbase), ddy(uvbase)), dayness ) * _Color;
				clip(texCol.a - _Cutoff);
				

				UNITY_LIGHT_ATTENUATION(attenuation, i, wPos);

				float3 specularTint;
				float oneMinusReflectivity;
				float smoothness = _Smoothness;
				float3 albedo = DiffuseAndSpecularFromMetallic(
					texCol, _Metallic, specularTint, oneMinusReflectivity
				);
				
				float3 viewDir = normalize(_WorldSpaceCameraPos - wPos);
				UnityLight light;
				light.color = attenuation * _LightColor0.rgb;
				light.dir = normalize(UnityWorldSpaceLightDir(i.wPos));
				UnityIndirect indirectLight;
				#ifdef UNITY_PASS_FORWARDADD
				indirectLight.diffuse = indirectLight.specular = 0;
				#else
				indirectLight.diffuse = max(0, ShadeSH9(float4(normal, 1)));
				float3 reflectionDir = reflect(-viewDir, normal);
				Unity_GlossyEnvironmentData envData;
				envData.roughness = 1 - smoothness;
				envData.reflUVW = reflectionDir;
				indirectLight.specular = Unity_GlossyEnvironment(
					UNITY_PASS_TEXCUBE(unity_SpecCube0), unity_SpecCube0_HDR, envData
				);
				#endif

				float3 col = UNITY_BRDF_PBS(
					albedo, specularTint,
					oneMinusReflectivity, smoothness,
					normal, viewDir,
					light, indirectLight
				);

				#ifdef UNITY_PASS_FORWARDADD
				return float4(col, 0);
				#else
				return float4(col, 1);
				#endif
			}
			#else
			float4 frag(v2f i) : SV_Target
			{
				float alpha = _Color.a;
				if (_Cutoff > 0)
					alpha *= tex2D(_MainTexDay, i.uv).a;
				clip(alpha - _Cutoff);
				SHADOW_CASTER_FRAGMENT(i)
			}
			#endif
		ENDCG

		Pass
		{
			Tags { "LightMode" = "ForwardBase" }
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma multi_compile_fwdbase_fullshadows
			#pragma multi_compile UNITY_PASS_FORWARDBASE
			ENDCG
		}

		Pass
		{
			Tags { "LightMode" = "ForwardAdd" }
			Blend One One
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma multi_compile_fwdadd_fullshadows
			#pragma multi_compile UNITY_PASS_FORWARDADD
			ENDCG
		}

		Pass
		{
			Tags { "LightMode" = "ShadowCaster" }
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma multi_compile_shadowcaster
			#pragma multi_compile UNITY_PASS_SHADOWCASTER
			ENDCG
		}
	}
}

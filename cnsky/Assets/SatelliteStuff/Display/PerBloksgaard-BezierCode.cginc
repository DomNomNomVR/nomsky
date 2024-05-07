// Copyright Per Bloksgaard, 2014 - https://perbloksgaard.dk
// I was inspired by https://www.shadertoy.com/view/XsX3zf but instead of a fast 
// distance approximation, I wanted the exact distance to a quadratic bezier spline.

//Blog article about this shader below. (In danish)
//http://www.hinnerup.net/permanent/2014/01/23/bezier_spline_shader/

#define dd(a) dot(a,a)
float addv(float2 a) { return a.x + a.y; }

//Find roots using Cardano's method. http://en.wikipedia.org/wiki/Cubic_function#Cardano.27s_method
float2 solveCubic2(float3 a)
{
	float p = a.y-a.x*a.x/3.;
	float p3 = p*p*p;
	float q = a.x*(2.*a.x*a.x-9.*a.y)/27.+a.z;
	float d = q*q+4.*p3/27.;
	//if( p < 10000 && d > 0) return 0;
	if(d>.0 )
	{
		// Original code from Per Bloksgaard has issues in a few cases here.
		float2 x = (float2(1,-1)*sqrt(d)-q)*.5;
		return (addv(sign(x)*pow(abs(x),(1./3.)))-a.x/3.).xx;
	}
	float v = acos(-sqrt(-27./p3)*q*.5)/3.;
	float m = cos(v);
	float n = sin(v)*sqrt(3);
	return float2(m+m,-n-m)*sqrt(-p/3.)-a.x/3.;
}

// How to solve the equation below can be seen on this image.
// http://www.perbloksgaard.dk/research/DistanceToQuadraticBezier.jpg
float calculateDistanceToQuadraticBezier(out float t, float2 p, float2 a, float2 b, float2 c)
{
	b += lerp((1e-4).xx,(0.).xx,abs(sign(b*2.-a-c)));
	float2 A = b-a;
	float2 B = c-b-A;
	float2 C = p-a;
	float2 D = A*2.;
	float2 T = clamp((solveCubic2(float3(-3.*dot(A,B),dot(C,B)-2.*dd(A),dot(C,A))/-dd(B))),0.,1.);
	
	// Added to know where you are along the line.
	float Ma = dd(C-(D+B*T.x)*T.x);
	float Mb = dd(C-(D+B*T.y)*T.y);
	t = (Ma<Mb)?T.x:T.y;
	return sqrt(min(Ma,Mb));
}

float calculateDistanceToQuadraticBezier3(out float t, float3 p, float3 a, float3 b, float3 c)
{
	b.xyz += lerp((1e-4).xxx,(0.).xxx,abs(sign(b.xyz*2.-a.xyz-c.xyz)));
	float3 A = b-a;
	float3 B = c-b-A;
	float3 C = p-a;
	float3 D = A*2.;
	
	
	// Figure out if we are gonna have a bad day.
	A.z *= 0.1;
	B.z *= 0.1;
	C.z *= 0.1;
	D.z *= 0.1;
	//GAHHHHHHHHHHHH PULL MY HAIR OUT
	float ddb = dd(B) + 0.0001;
	float2 T = clamp((solveCubic2(float3(-3.*dot(A,B),dot(C,B)-2.*dd(A),dot(C,A))/-ddb)),0.,1.);

	// Added to know where you are along the line.
	float Ma = dd(C-(D+B*T.x)*T.x);
	float Mb = dd(C-(D+B*T.y)*T.y);
	t = (Ma<Mb)?T.x:T.y;

	return sqrt(min(Ma,Mb));
}

float cross2d( float2 A, float2 B )
{
	return A.x * B.y - A.y * B.x;
}

void ResolveBezierGeometry( out float3 viewspacepos[5], float3 bez[3], float expand)
{
	float3 base[3] = {
		(bez[0]),
		(bez[1]),
		(bez[2]) };

	float3 clip2d[3] = {
		base[0].xyz,
		base[1].xyz,
		base[2].xyz,
	};

	// Axis to expand from: In X-Y, 
	
	float3 orthos[3];
	orthos[0] = normalize( cross( clip2d[0] - clip2d[1], clip2d[0] ) ); //(clip2d[0] - clip2d[1]).yx * float2( 1.0, -1.0 );
	orthos[2] = normalize( cross( clip2d[2] - clip2d[1], clip2d[2] ) ); //(clip2d[2] - clip2d[1]).yx * float2( 1.0, -1.0 );
	orthos[0] = normalize( orthos[0] )*expand;
	orthos[2] = normalize( orthos[2] )*expand;
	
	// TODO: Try to expand from point on line from A--B orthogonal to center point.
	orthos[1] = normalize( orthos[0] - orthos[2] ) * expand;
		
	float clockwisiness = sign( cross2d( clip2d[2].xy - clip2d[1].xy, clip2d[0].xy - clip2d[1].xy ) );
	
	viewspacepos[0] = base[0] - float3( orthos[0] ) * clockwisiness;
	viewspacepos[1] = base[0] + float3( orthos[0] ) * clockwisiness;
	viewspacepos[2] = base[1] - float3( orthos[1] ) * clockwisiness;
	viewspacepos[3] = base[2] - float3( orthos[2] ) * clockwisiness;
	viewspacepos[4] = base[2] + float3( orthos[2] ) * clockwisiness;
}

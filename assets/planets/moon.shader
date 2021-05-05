shader_type canvas_item;

uniform vec4 planet_color : hint_color = vec4(0.78, 0.85, 1.0, 1.0);
uniform vec4 atmosphere_color : hint_color = vec4(0.93, 0.93, 1.0, 1.0);
uniform float atmosphere_strength : hint_range(0, 1) = 1.0;
uniform float planet_size : hint_range(0, 1) = 0.81;
uniform float planet_ruggedness : hint_range(0, 1) = 0.05;
uniform float rotation_speed : hint_range(0, 1) = 0.033;

uniform vec3 light_direction = vec3(0.7, -0.5, 0.5);
uniform float light_wrap : hint_range(0, 1) = 0.33;

uniform float in_inner = 0.2;
uniform float in_outer = 0.2;
uniform float out_inner = 0.2;
uniform float out_outer = 0.5;

uniform float antialiasing : hint_range(1, 5, 1) = 2.0;

float noise3D(vec3 p) {
	return fract(sin(dot(p, vec3(12.9898,78.233,128.852))) * 43758.5453)*2.0-1.0;
}

float simplex3D(vec3 p) {
	float f3 = 1.0/3.0;
	float s = (p.x+p.y+p.z)*f3;
	int i = int(floor(p.x+s));
	int j = int(floor(p.y+s));
	int k = int(floor(p.z+s));

	float g3 = 1.0/6.0;
	float t = float((i+j+k))*g3;
	float x0 = float(i)-t;
	float y0 = float(j)-t;
	float z0 = float(k)-t;
	x0 = p.x-x0;
	y0 = p.y-y0;
	z0 = p.z-z0;
	int i1,j1,k1;
	int i2,j2,k2;
	if(x0>=y0) {
		if		(y0>=z0){ i1=1; j1=0; k1=0; i2=1; j2=1; k2=0; } // X Y Z order
		else if	(x0>=z0){ i1=1; j1=0; k1=0; i2=1; j2=0; k2=1; } // X Z Y order
		else 			{ i1=0; j1=0; k1=1; i2=1; j2=0; k2=1; } // Z X Z order
	} else {
		if		(y0<z0) { i1=0; j1=0; k1=1; i2=0; j2=1; k2=1; } // Z Y X order
		else if	(x0<z0) { i1=0; j1=1; k1=0; i2=0; j2=1; k2=1; } // Y Z X order
		else 			{ i1=0; j1=1; k1=0; i2=1; j2=1; k2=0; } // Y X Z order
	}
	float x1 = x0 - float(i1) + g3;
	float y1 = y0 - float(j1) + g3;
	float z1 = z0 - float(k1) + g3;
	float x2 = x0 - float(i2) + 2.0*g3;
	float y2 = y0 - float(j2) + 2.0*g3;
	float z2 = z0 - float(k2) + 2.0*g3;
	float x3 = x0 - 1.0 + 3.0*g3;
	float y3 = y0 - 1.0 + 3.0*g3;
	float z3 = z0 - 1.0 + 3.0*g3;
	vec3 ijk0 = vec3(float(i), float(j), float(k));
	vec3 ijk1 = vec3(float(i+i1), float(j+j1), float(k+k1));
	vec3 ijk2 = vec3(float(i+i2), float(j+j2), float(k+k2));
	vec3 ijk3 = vec3(float(i+1), float(j+1), float(k+1));
	vec3 gr0 = normalize(vec3(noise3D(ijk0),noise3D(ijk0*2.01),noise3D(ijk0*2.02)));
	vec3 gr1 = normalize(vec3(noise3D(ijk1),noise3D(ijk1*2.01),noise3D(ijk1*2.02)));
	vec3 gr2 = normalize(vec3(noise3D(ijk2),noise3D(ijk2*2.01),noise3D(ijk2*2.02)));
	vec3 gr3 = normalize(vec3(noise3D(ijk3),noise3D(ijk3*2.01),noise3D(ijk3*2.02)));
	float n0 = 0.0;
	float n1 = 0.0;
	float n2 = 0.0;
	float n3 = 0.0;
	float t0 = 0.5 - x0*x0 - y0*y0 - z0*z0;
	if(t0>=0.0) {
		t0*=t0;
		n0 = t0 * t0 * dot(gr0, vec3(x0, y0, z0));
	}
	float t1 = 0.5 - x1*x1 - y1*y1 - z1*z1;
	if(t1>=0.0) {
		t1*=t1;
		n1 = t1 * t1 * dot(gr1, vec3(x1, y1, z1));
	}
	float t2 = 0.5 - x2*x2 - y2*y2 - z2*z2;
	if(t2>=0.0) {
		t2 *= t2;
		n2 = t2 * t2 * dot(gr2, vec3(x2, y2, z2));
	}
	float t3 = 0.5 - x3*x3 - y3*y3 - z3*z3;
	if(t3>=0.0) {
		t3 *= t3;
		n3 = t3 * t3 * dot(gr3, vec3(x3, y3, z3));
	}
	return 96.0*(n0+n1+n2+n3);
}

float fbm(vec3 p) {
	float f;
	f  = 0.50000*(simplex3D(p)); p = p*2.01;
	f += 0.25000*(simplex3D(p)); p = p*2.02;
	f += 0.12500*(simplex3D(p)); p = p*2.03;
	f += 0.06250*(simplex3D(p)); p = p*2.04;
	f += 0.03125*(simplex3D(p)); p = p*2.05;
	f += 0.015625*(simplex3D(p));
	return f;
}

vec3 rotate(vec3 p, float theta) {
	return mat3(
			vec3(cos(theta), 0, sin(theta)),
			vec3(0, 1, 0),
			vec3(-sin(theta), 0, cos(theta))) * p;
}

vec4 pixel(vec2 uv, float time) {
	// planet centered on 0.5, 0.5, so now pos is -0.5 to 0.5
	vec3 pos = vec3(uv.xy - 0.5, 0.0);

	// LIGHT
	//vec3 l = normalize(vec3(sin(TIME*0.1), cos(TIME*0.1), 0));
	vec3 light = normalize(light_direction);

	// PLANET
	float rad = planet_size / 2.0; // radius
	float z_in = sqrt(rad*rad - pos.x*pos.x - pos.y*pos.y);
	float z_out = sqrt(-rad*rad + pos.x*pos.x + pos.y*pos.y);

	vec3 sphere_pos = normalize(vec3(pos.x, pos.y, z_in));
	//return vec4(sphere_pos, 1);

	// NORMALS
	vec3 norm = normalize(vec3(pos.x, pos.y, z_in)); // normals from sphere
	norm = rotate(norm, time * rotation_speed);
	vec3 norm_out = normalize(vec3(pos.x, pos.y, z_out)); // normals from outside sphere
	float e = planet_ruggedness; // planet rugosity
	float nx = fbm(vec3(norm.x+e, norm.y,   norm.z  ))*0.5+0.5; // x normal displacement
	float ny = fbm(vec3(norm.x,   norm.y+e, norm.z  ))*0.5+0.5; // y normal displacement
	float nz = fbm(vec3(norm.x,   norm.y,   norm.z+e))*0.5+0.5; // z normal displacement
	norm = normalize(vec3(norm.x*nx, norm.y*ny, norm.z*nz));
	//norm = (norm+1.)/2.; // for normals visualization
	//return vec4(norm, 1);
	
	// TEXTURE
	float n = 1.0-(fbm(norm)*0.5+0.5); // noise for every pixel in planet
	//return vec4(n,n,n,1);
	// ATMOS
	float z_in_atm  = (rad * in_outer)  / z_in - in_inner;   // inner atmos
	float z_out_atm = (rad * out_inner) / z_out - out_outer; // outer atmos
	z_in_atm = max(0.0, z_in_atm);
	z_out_atm = max(0.0, z_out_atm);

	// DIFFUSE LIGHT
	float diffuse = max(0.0, dot(sphere_pos, light));
	float diffuse_out = max(0.0, dot(norm_out, light) + light_wrap);
	//return vec4(vec3(diffuse_out), 1);

	//return vec4(vec3(n * diffuse),1.0);
	//return vec4(vec3(z_in_atm * diffuse),1.0);
	//return vec4(vec3(z_out_atm * diffuse_out),1.0);
	
	float atm = z_in_atm * diffuse + z_out_atm * diffuse_out;
	
	return vec4(planet_color.rgb * n * diffuse, 1.0) + vec4(atmosphere_color.rgb * atm, atmosphere_strength);
}

void fragment() {
	vec2 o = TEXTURE_PIXEL_SIZE * (0.125 / antialiasing);
	vec4 color = vec4(0.0);
	for(float x = -antialiasing; x <= antialiasing + 0.01; x++) {
		for(float y = -antialiasing; y <= antialiasing + 0.01; y++) {
			color += pixel(UV + vec2(o.x * x, o.y * y), TIME);
		}
	}
	COLOR = color / ((antialiasing * 2.0 + 1.0) * (antialiasing * 2.0 + 1.0));
	//COLOR = pixel(UV, TIME);
}
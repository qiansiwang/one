/*
One is designed to programatically create 3D points.
The end result is one dimentional array buffer.
*/

/*
Origin is the base class any 3D object
x, y, z translates, i, j, k rotates
either pass [s, i, j, k] as quaternion to rotation
or pass [theta, i, j, k] as an angle and a vector
*/
class Origin {
	constructor (translation, rotation, option){
		this.translation = translation;
		this.rotation = rotation;
		this.option = option;
	};
	get s(){
		return this.rotation[0];
	}
	get theta(){
		return this.rotation[0];
	}
	get i(){
		return this.rotation[1];
	}
	get j(){
		return this.rotation[2];
	}
	get k(){
		return this.rotation[3];
	}
	rotateByQuaternion(){
		let e11 = 2*(this.s*this.s + this.i*this.i) - 1;
		let e12 = 2*(this.i*this.j - this.s*this.k);
		let e13 = 2*(this.i*this.k + this.s*this.j);
		let e21 = 2*(this.i*this.j + this.s*this.k);
		let e22 = 2*(this.s*this.s + this.j*this.j) - 1;
		let e23 = 2*(this.j*this.k - this.s*this.i);
		let e31 = 2*(this.i*this.k - this.s*this.j);
		let e32 = 2*(this.j*this.k + this.s*this.i);
		let e33 = 2*(this.s*this.s + this.k*this.k) - 1;
		let x = e11*this.translation[0] + e12*this.translation[1] + e13*this.translation[2];
		let y = e21*this.translation[0] + e22*this.translation[1] + e23*this.translation[2];
		let z = e31*this.translation[0] + e32*this.translation[1] + e33*this.translation[2];
		return [x, y, z];
	}
	rotateByAnglarVector(){
		let cos = Math.cos(2 * this.theta);
		let sin = Math.sin(2 * this.theta);
		let vx = this.normalize(this.translation[0], this.translation[1], this.translation[2])[0];
		let vy = this.normalize(this.translation[0], this.translation[1], this.translation[2])[1];
		let vz = this.normalize(this.translation[0], this.translation[1], this.translation[2])[2];
		let ux = this.normalize(this.i, this.j, this.k)[0];
		let uy = this.normalize(this.i, this.j, this.k)[1];
		let uz = this.normalize(this.i, this.j, this.k)[2];
		let vudot = vx*ux + vy*uy + vz*uz;
		let newx = (1-cos)*vudot*vx + cos*ux + sin*(vy*uz-vz*uy);
		let newy = (1-cos)*vudot*vy + cos*uy + sin*(-vx*uz+vz*ux*vz);
		let newz = (1-cos)*vudot*vz	+ cos*uz + sin*(vx*uy-vy*ux);
		return [newx, newy, newz];
	}
	get X(){
		if (this.option == "A"){
			return this.rotateByAnglarVector[0];	
		}
		else {
			return this.rotateByQuaternion[0];
		}
	}
	get Y(){
		if (this.option == "A"){
			return this.rotateByAnglarVector[1];	
		}
		else {
			return this.rotateByQuaternion[1];
		}
	}
	get Z(){
		if (this.option == "A"){
			return this.rotateByAnglarVector[2];
		}
		else {
			return this.rotateByQuaternion[2];
		}
	}		
	static normalize(x, y, z){
		let a = Math.sqrt(x*x + y*y + z*z)
		//0,0,0 will give infinity
		return [x/a, y/a, z/a]
	}
};

/*
Single Point is single vertex with origin
*/

class SinglePoint extends Origin {
	constructor (vector, translation, rotation, option) {
		super([vector[0]+translation[0], vector[1]+translation[1], vector[2]+translation[2]], rotation, option);
		this.localVector = vector
	}
    get arrayBuffer(){
		return [this.X, this.Y, this.Z];
	}
	get originDistance() {
		return Math.sqrt(this.X*this.X+this.Y*this.Y+this.Z*this.Z);
	}
}


/*
Segment is consisted of two points. If two points are same return undefined
array buffer return A, B
*/

class Segment extends Origin {
	constructor (vectora, vectorb, translation, rotation, option){
		if (vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]){
			return;
		}
		else {
			super(translation, rotation, option);
			this.vectora = vectora;
			this.vectorb = vectorb;
		}
	}
	get pointA(){
		return new SinglePoint(this.vectora, this.translation, this.rotation, this.option);
	}
	get pointB(){
		return new SinglePoint(this.vectorb, this.translation, this.rotation, this.option);
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z];
	}
    get segmentVector(){
		var x = this.pointB.X - this.pointA.X;
		var y = this.pointB.Y - this.pointA.Y;
		var z = this.pointB.Z - this.pointA.Z;
		return [x, y, z];
	}
	get vectorA(){
		return this.pointA.arrayBuffer;
	}
	get vectorB(){
		return this.pointB.arrayBuffer;
	}
	collinear (x, y, z){
		let x0 = this.pointA.X;
		let y0 = this.pointA.Y;
		let z0 = this.pointA.Z;
		let a = this.segmentVector[0];
		let b = this.segmentVector[1];
		let c = this.segmentVector[2];
		let t1, t2, t3
		t1 = (x-x0)/a;
		t2 = (y-y0)/b;
		t3 = (z-z0)/c;
		if (t1 === t2 && t2 === t3){
			return true;
		}
		else {
			return false;
		}
	}
}

/*
Tria is triangular shape. If any two points are same return undefined
If three points are collinear return undefined
*/

class Tria extends Origin {
	constructor (vectora, vectorb, vectorc, translation, rotation, option){
		if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2])){
			return;
		}
		let tempseg = new Segment(vectora, vectorb, [0,0,0], [0,0,0,0], "A")
		if (tempseg.collinear(vectorc[0], vectorc[1], vectorc[2])){
			return;
		}
		else {
			super(translation, rotation, option);
			this.vectora = vectora;
			this.vectorb = vectorb;
			this.vectorc = vectorc;
		}
	}
	get pointA(){
		return new SinglePoint(this.vectora, this.translation, this.rotation, this.option);
	}
	get pointB(){
		return new SinglePoint(this.vectorb, this.translation, this.rotation, this.option);
	}
	get pointC(){
		return new SinglePoint(this.vectorc, this.translation, this.rotation, this.option);
	}
	get segmentBC(){
		return new Segment(this.vectorb, this.vectorc, this.translation, this.rotation, this.option);
	}
	get segmentAB(){
		return new Segment(this.vectora, this.vectorb, this.translation, this.rotation, this.option);
	}
	get segmentCA(){
		return new Segment(this.vectorc, this.vectora, this.translation, this.rotation, this.option);
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z, this.pointC.X, this.pointC.Y, this.pointC.Z];
	}
	get perpendicularVector(){
		let ab = this.segmentAB.segmentVector;
		let bc = this.segmentBC.segmentVector;
		let i = ab[1]*bc[2] - ab[2]*bc[1];
		let j = ab[2]*bc[0] - ab[0]*bc[2];
		let k = ab[0]*bc[1] - ab[1]*bc[0];
		return [i, j, k];
	}

	coplanar(x, y, z){
		let a = this.perpendicularVector[0];
		let b = this.perpendicularVector[1];
		let c = this.perpendicularVector[2];
		let x0 = this.pointA.X;
		let y0 = this.pointA.Y;
		let z0 = this.pointA.Z;
		let result = a*(x-x0) + b*(y-y0) + c*(z-z0);
		if (result === 0) {
			return true;
		}
		else {
			return false;
		}
	}
}

/*
Quad is four points on same surface. If any two points are same return undefined.
If the four points are not coplanar then return undefined.
*/

class Quad extends Tria {
	constructor (vectora, vectorb, vectorc, vectord, translation, rotation, option){
		if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2]) ||
			(vectorc[0] === vectord[0] && vectorc[1] === vectord[1] && vectorc[2] === vectord[2]) ||
			(vectora[0] === vectord[0] && vectora[1] === vectord[1] && vectora[2] === vectord[2]) ||
			(vectorb[0] === vectord[0] && vectorb[1] === vectord[1] && vectorb[2] === vectord[2])
			){
			return;
		}
		let temptri = new Tria(vectora, vectorb, vectorc, [0,0,0], [0,0,0,0], "A");
		if (temptri.coplanar(vectord[0], vectord[1], vectord[2])){
			super(vectora, vectorb, vectorc, translation, rotation, option);
			this.vectord = vectord;
		}
		else {
			return;
		}
	}
	get pointD(){
		return new SinglePoint(this.vectord, this.translation, this.rotation, this.option);
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z, this.pointC.X, this.pointC.Y, this.pointC.Z, 
		this.pointD.X, this.pointD.Y, this.pointD.Z];
	}
	get segmentCD(){
		return new Segment(this.vectorc, this.vectord, this.translation, this.rotation, this.option);
	}
	get segmentDA(){
		return new Segment(this.vectord, this.vectora, this.translation, this.rotation, this.option);
	}
}

/*
Orthogonal Box is consisted of two Quads
Origin stays in center 
*/

class Box extends Origin {
	constructor (length, width, height, translation, rotation, option){
		super(translation, rotation, option);
		this.length = length;
		this.width = width;
		this.height = height;
	}
	get topQuad(){
		let A = [- this.width/2, this.height/2, - this.length/2];
		let B = [this.width/2, this.height/2, - this.length/2];
		let C = [this.width/2, this.height/2, this.length/2];
		let D = [- this.width/2, this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation, this.option);
		return result;
	}
	get bottomQuad(){
		let A = [- this.width/2, - this.height/2, - this.length/2];
		let B = [this.width/2, - this.height/2, - this.length/2];
		let C = [this.width/2, - this.height/2, this.length/2];
		let D = [- this.width/2, - this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation, this.option);
		return result;
	}
	get arrayBuffer(){
		let result = []
		this.topQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		this.bottomQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		return result;
	}
}

/*
Polygons are made of lots of segment
Polygons have one center and translation vector from origin point
If origin is 0, 0, 0 assuming translation vector j
Polygons plane is perpendicular to translation vector
*/

class Polygon extends Origin {
	constructor (edge, radius, centervector, translation, rotation, option){
		super(translation, rotation, option);
		this.edge = edge;
		this.radius = radius;
        this.center = centervector;
	}
	get centerPoint(){
		return new SinglePoint(this.center, this.translation, this.rotation, this.option);
	}
	
	
	
	
}

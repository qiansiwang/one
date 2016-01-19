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
	constructor (translation, rotation){
		this.translation = translation;
		this.rotation = rotation;
	};
	get option(){
		return this.rotation[4];
	}
	get s(){
		if (this.option === "A") {
			return;
		}
		else {
			return this.rotation[0];
		}
	}
	get theta(){
		if (this.option === "A") {
			return this.rotation[0]/180*Math.PI/2;
		}
		else {
			return;
		}
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
	get rotateByQuaternion(){
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
	get rotateByAnglarVector(){
		let cos = Math.cos(2 * this.theta);
		let sin = Math.sin(2 * this.theta);
		let ux = this.translation[0];
		let uy = this.translation[1];
		let uz = this.translation[2];
		let normal2 = this.normalize(this.i, this.j, this.k)
		let vx = normal2[0];
		let vy = normal2[1];
		let vz = normal2[2];
		let vudot = vx*ux + vy*uy + vz*uz;
		let newx = (1-cos)*vudot*vx + cos*ux + sin*(vy*uz-vz*uy);
		let newy = (1-cos)*vudot*vy + cos*uy + sin*(-vx*uz+vz*ux*vz);
		let newz = (1-cos)*vudot*vz	+ cos*uz + sin*(vx*uy-vy*ux);
		return [newx, newy, newz];
	}
	get X(){
		if (this.option === "A"){
			return this.rotateByAnglarVector[0];
		}
		else {
			return this.rotateByQuaternion[0];
		}
	}
	get Y(){
		if (this.option === "A"){
			return this.rotateByAnglarVector[1];
		}
		else {
			return this.rotateByQuaternion[1];
		}
	}
	get Z(){
		if (this.option === "A"){
			return this.rotateByAnglarVector[2];
		}
		else {
			return this.rotateByQuaternion[2];
		}
	}
	normalize(x, y, z){
		let a = Math.sqrt(x*x + y*y + z*z)
		//0,0,0 will give infinity
		if (x==y && y==z && z==0) {
			throw "Rotation Vector Cannot Be (0,0,0)"
		}
		return [x/a, y/a, z/a]
	}
};

/*
Single Point is single vertex with origin
*/

class SinglePoint extends Origin {
	constructor (vector, translation, rotation) {
		super([vector[0]+translation[0], vector[1]+translation[1], vector[2]+translation[2]], rotation);
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
	constructor (vectora, vectorb, translation, rotation){
		if (vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]){
			return;
		}
		else {
			super(translation, rotation);
			this.vectora = vectora;
			this.vectorb = vectorb;
		}
	}
	get pointA(){
		return new SinglePoint(this.vectora, this.translation, this.rotation);
	}
	get pointB(){
		return new SinglePoint(this.vectorb, this.translation, this.rotation);
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z];
	}
    get segmentVector(){
		let pointa = this.pointA
		let pointb = this.pointB
		var x = pointb.X - pointa.X;
		var y = pointb.Y - pointa.Y;
		var z = pointb.Z - pointa.Z;
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
	constructor (vectora, vectorb, vectorc, translation, rotation){
		if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2])){
			return;
		}
		let tempseg = new Segment(vectora, vectorb, [0,0,0], [0,1,0,0,"A"])
		if (tempseg.collinear(vectorc[0], vectorc[1], vectorc[2])){
			return;
		}
		else {
			super(translation, rotation);
			this.vectora = vectora;
			this.vectorb = vectorb;
			this.vectorc = vectorc;
		}
	}
	get pointA(){
		return new SinglePoint(this.vectora, this.translation, this.rotation);
	}
	get pointB(){
		return new SinglePoint(this.vectorb, this.translation, this.rotation);
	}
	get pointC(){
		return new SinglePoint(this.vectorc, this.translation, this.rotation);
	}
	get segmentBC(){
		return new Segment(this.vectorb, this.vectorc, this.translation, this.rotation);
	}
	get segmentAB(){
		return new Segment(this.vectora, this.vectorb, this.translation, this.rotation);
	}
	get segmentCA(){
		return new Segment(this.vectorc, this.vectora, this.translation, this.rotation);
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z, this.pointC.X, this.pointC.Y, this.pointC.Z];
	}
	get perpendicularVector(){
		let segab = this.segmentAB;
		let ab = segab.segmentVector
		let segbc = this.segmentBC;
		let bc = segbc.segmentVector
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

class Quad extends Origin {
	constructor (vectora, vectorb, vectorc, vectord, translation, rotation){
		if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2]) ||
			(vectorc[0] === vectord[0] && vectorc[1] === vectord[1] && vectorc[2] === vectord[2]) ||
			(vectora[0] === vectord[0] && vectora[1] === vectord[1] && vectora[2] === vectord[2]) ||
			(vectorb[0] === vectord[0] && vectorb[1] === vectord[1] && vectorb[2] === vectord[2])
			){
			return;
		}
		let temptri = new Tria(vectora, vectorb, vectorc, [0,0,0], [0,1,0,0,"A"]);
		if (temptri.coplanar(vectord[0], vectord[1], vectord[2])){
			super(translation, rotation);
			this.vectora = vectora;
			this.vectorb = vectorb;
			this.vectorc = vectorc;
			this.vectord = vectord;
		}
		else {
			return;
		}
	}
	get triaA(){
		return new Tria(this.vectora, this.vectorb, this.vectorc, this.translation, this.rotation)
	}
	get triaB(){
		return new Tria(this.vectorc, this.vectord, this.vectora, this.translation, this.rotation)
	}
	get pointA(){
		return new SinglePoint(this.vectora, this.translation, this.rotation);
	}
	get pointB(){
		return new SinglePoint(this.vectorb, this.translation, this.rotation);
	}
	get pointC(){
		return new SinglePoint(this.vectorc, this.translation, this.rotation);
	}
	get pointD(){
		return new SinglePoint(this.vectord, this.translation, this.rotation);
	}
	get arrayBuffer(){
		return this.triaA.arrayBuffer.concat(this.triaB.arrayBuffer)
	}
	get segmentCD(){
		return new Segment(this.vectorc, this.vectord, this.translation, this.rotation);
	}
	get segmentDA(){
		return new Segment(this.vectord, this.vectora, this.translation, this.rotation);
	}
}

/*
Orthogonal Box is consisted of two Quads
Origin stays in center
*/

class Box extends Origin {
	constructor (length, width, height, translation, rotation){
		super(translation, rotation);
		this.length = length;
		this.width = width;
		this.height = height;
	}
	get topQuad(){
		let A = [- this.width/2, this.height/2, - this.length/2];
		let B = [- this.width/2, this.height/2, this.length/2];
		let C = [this.width/2, this.height/2, this.length/2];
		let D = [this.width/2, this.height/2, - this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation);
		return result;
	}
	get leftQuad(){
		let A = [- this.width/2, this.height/2, - this.length/2];
		let B = [- this.width/2, - this.height/2, - this.length/2];
		let C = [- this.width/2, - this.height/2, this.length/2];
		let D = [- this.width/2, this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation);
		return result;
	}
	get frontQuad(){
		let A = [- this.width/2, this.height/2, this.length/2];
		let B = [- this.width/2, - this.height/2, this.length/2];
		let C = [this.width/2, - this.height/2, this.length/2];
		let D = [this.width/2, this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation);
		return result;
	}
	get bottomQuad(){
		let A = [- this.width/2, - this.height/2, - this.length/2];
		let B = [this.width/2, - this.height/2, - this.length/2];
		let C = [this.width/2, - this.height/2, this.length/2];
		let D = [- this.width/2, - this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation);
		return result;
	}
	get rightQuad(){
		let A = [this.width/2, this.height/2, this.length/2];
		let B = [this.width/2, - this.height/2, this.length/2];
		let C = [this.width/2, - this.height/2, - this.length/2];
		let D = [this.width/2, this.height/2, - this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation);
		return result;
	}
	get backQuad(){
		let A = [this.width/2, this.height/2, - this.length/2];
		let D = [- this.width/2, this.height/2, - this.length/2];
		let C = [- this.width/2, - this.height/2, - this.length/2];
		let B = [this.width/2, - this.height/2, - this.length/2];
		let result = new Quad(A,B,C,D,this.translation, this.rotation);
		return result;
	}
	get arrayBuffer(){
		let result = []
		this.topQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		this.leftQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		this.frontQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		this.bottomQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		this.rightQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		this.backQuad.arrayBuffer.forEach(function (e){
			result.push(e);
		})
		return result;
	}
}

/*
Polygons are made of lots of triangulars parallel to XZ plane
Polygons have one center and translation vector from origin point
*/

class Polygon extends Origin {
	constructor (edge, radius, centervector, translation, rotation){
		super(translation, rotation);
		this.edge = edge;
		this.radius = radius;
    this.center = centervector;
	}
	get centerPoint(){
		return new SinglePoint(this.center, this.translation, this.rotation);
	}
	get triaList(){
		let result = [];
		for (let i = 0;	i	<	this.edge; i++){
			let center = this.center;
			let a = Math.PI*2/this.edge;
			let x1 = Math.sin(a*i)*this.radius + this.center[0];
			let y1 = this.center[1];
			let z1 = Math.cos(a*i)*this.radius + this.center[2];
			let point1 = [x1, y1, z1];
			let x2 = Math.sin(a*(i+1))*this.radius + this.center[0];
			let y2 = this.center[1];
			let z2 = Math.cos(a*(i+1))*this.radius + this.center[2];
			let point2 = [x2, y2, z2];
			let tria = new Tria(center, point1, point2, this.translation, this.rotation);
			result.push(tria);
		}
		return result;
	}
	get arrayBuffer(){
		let result = [];
		let list = this.triaList;
		list.forEach(function(e){
			let x0 = e.pointA.X;
			let y0 = e.pointA.Y;
			let z0 = e.pointA.Z;
			let x1 = e.pointB.X;
			let y1 = e.pointB.Y;
			let z1 = e.pointB.Z;
			let x2 = e.pointC.X;
			let y2 = e.pointC.Y;
			let z2 = e.pointC.Z
			result.push(x0);
			result.push(y0);
			result.push(z0);
			result.push(x1);
			result.push(y1);
			result.push(z1);
			result.push(x2);
			result.push(y2);
			result.push(z2);
		})
		return result;
	}
}

/*
Orthogonal Cylinder
*/
class Cylinder extends Origin {
	constructor (edge, radius, height, translation, rotation){
		super(translation, rotation);
		this.edge = edge;
		this.radius = radius;
		this.height = height;
	}
	get topPolygon(){
		return new Polygon(this.edge, this.radius, [0,this.height/2,0], this.translation, this.rotation);
	}
	get bottomPolygon(){
		return new Polygon(this.edge, this.radius, [0,-this.height/2,0], this.translation, this.rotation);
	}
	get arrayBuffer(){
		return this.topPolygon.arrayBuffer.concat(this.bottomPolygon.arrayBuffer);
	}
}

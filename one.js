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
		let newy = (1-cos)*vudot*vy + cos*uy + sin*(-vx*uz+vz*ux);
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
		super([vector[0], vector[1], vector[2]], rotation);
		this.localVector = vector
	}
    get arrayBuffer(){
		let result = new Float32Array(3)
		result[0] = this.X + this.translation[0];
		result[1] = this.Y + this.translation[1];
		result[2] = this.Z + this.translation[2];
		return result;
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
		/*if (vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]){
			return;
		}
		else {
		}*/
		super(translation, rotation);
		this.vectora = vectora;
		this.vectorb = vectorb;

	}
	get pointA(){
		return new SinglePoint(this.vectora, this.translation, this.rotation);
	}
	get pointB(){
		return new SinglePoint(this.vectorb, this.translation, this.rotation);
	}
	get arrayBuffer(){
		let result = new Float32Array(6)
		result[0] = this.pointA.X + this.translation[0]
		result[1] = this.pointA.Y + this.translation[1]
		result[2] = this.pointA.Z + this.translation[2]
		result[3] = this.pointB.X + this.translation[0]
		result[4] = this.pointB.Y + this.translation[1]
		result[5] = this.pointB.Z + this.translation[2]
		return result;
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
		/*if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2])){
			return;
		}
		let tempseg = new Segment(vectora, vectorb, [0,0,0], [0,1,0,0,"A"])
		if (tempseg.collinear(vectorc[0], vectorc[1], vectorc[2])){
			return;
		}
		else {
		}*/
		super(translation, rotation);
		this.vectora = vectora;
		this.vectorb = vectorb;
		this.vectorc = vectorc;
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
		let result = new Float32Array(9);
		result[0] = this.pointA.X + this.translation[0]
		result[1] = this.pointA.Y + this.translation[1]
		result[2] = this.pointA.Z + this.translation[2]
		result[3] = this.pointB.X + this.translation[0]
		result[4] = this.pointB.Y + this.translation[1]
		result[5] = this.pointB.Z + this.translation[2]
		result[6] = this.pointC.X + this.translation[0]
		result[7] = this.pointC.Y + this.translation[1]
		result[8] = this.pointC.Z + this.translation[2]
		return result;
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
Not checking this any more. If the four points are not coplanar then return undefined.
*/

class Quad extends Origin {
	constructor (vectora, vectorb, vectorc, vectord, translation, rotation){
		/*if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
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

		}
		else {
			return;
		}*/
		super(translation, rotation);
		this.vectora = vectora;
		this.vectorb = vectorb;
		this.vectorc = vectorc;
		this.vectord = vectord;
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
		let a = this.triaA.arrayBuffer
		let b = this.triaB.arrayBuffer
		let result = new Float32Array(18);
		for (let i = 0; i < 18; i++){
			if (i < 9){
				result[i] = a[i];
			}
			else {
				result[i] = b[i-9];
			}
		}
		return result;
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
		let result = new Float32Array(108);
		for (let i=0; i<108; i++){
			let j = i%18;
			if (i < 18){
				result[i] = this.topQuad.arrayBuffer[j];
			}
			else if (i>= 18 && i<36){
				result[i] = this.leftQuad.arrayBuffer[j];
			}
			else if (i>=36 && i<54){
				result[i] = this.frontQuad.arrayBuffer[j];
			}
			else if (i>=54 && i<72){
				result[i] = this.bottomQuad.arrayBuffer[j];
			}
			else if (i>=72 && i<90){
				result[i] = this.rightQuad.arrayBuffer[j];
			}
			else {
				result[i] = this.backQuad.arrayBuffer[j];
			}
		}
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
		for (let i = 0;	i<this.edge; i++){
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
		let result = new Float32Array(this.edge*9);
		let list = this.triaList;
		for (let i=0; i< list.length; i++){
			let e = list[i];
			let x0 = e.pointA.X;
			let y0 = e.pointA.Y;
			let z0 = e.pointA.Z;
			let x1 = e.pointB.X;
			let y1 = e.pointB.Y;
			let z1 = e.pointB.Z;
			let x2 = e.pointC.X;
			let y2 = e.pointC.Y;
			let z2 = e.pointC.Z;
			result[i*9] = x0;
			result[i*9+1] = y0;
			result[i*9+2] = z0;
			result[i*9+3] = x1;
			result[i*9+4] = y1;
			result[i*9+5] = z1;
			result[i*9+6] = x2;
			result[i*9+7] = y2;
			result[i*9+8] = z2;
		}
		return result;
	}
}

/*
Cylinder takes top and bottom radius, and height.
*/
class Cylinder extends Origin {
	constructor (edge, topradius, bottomradius, height, translation, rotation){
		super(translation, rotation);
		this.edge = edge;
		this.topradius = topradius;
		this.bottomradius = bottomradius;
		this.height = height;
	}
	get topPolygon(){
		return new Polygon(this.edge, this.topradius, [0,this.height/2,0], this.translation, this.rotation);
	}
	get bottomPolygon(){
		return new Polygon(this.edge, this.bottomradius, [0,-this.height/2,0], this.translation, this.rotation);
	}
	get cylindricalSurface(){
		let result = new Float32Array(18*this.edge);
		for (let i = 0; i<this.edge; i++){
			let triatop = this.topPolygon.triaList[i];
			let triabottom = this.bottomPolygon.triaList[i];
			let a = triatop.vectorb;
			let b = triabottom.vectorb;
			let c = triabottom.vectorc;
			let d = triatop.vectorc;
			let quad = new Quad(a,b,c,d,this.translation,this.rotation);
			result.set(quad.arrayBuffer,i*18)
		}
		return result
	}
	get arrayBuffer(){
		let topbuffer = this.topPolygon.arrayBuffer;
		let bottombuffer = this.bottomBuffer;
		let cyli = this.cylindricalSurface;
		let result = new Float32Array(topbuffer.length + bottombuffer.length + cyli.length);
		result.set(topbuffer,0);
		result.set(cyli,topbuffer.length);
		result.set(bottombuffer,topbuffer.length + cyli.length);
		return result;
	}
	get bottomBuffer(){
		let bottombuffer = this.bottomPolygon.arrayBuffer;
		let bottomlength = bottombuffer.length
		let result = new Float32Array(bottombuffer.length);
		for (let i= bottomlength/3 -1; i>= 0; i--){
			let x = bottombuffer[i*3];
			let y = bottombuffer[i*3 + 1];
			let z = bottombuffer[i*3 + 2];
			result[bottomlength-i*3-3] = x;
			result[bottomlength-i*3-2] = y;
			result[bottomlength-i*3-1] = z;
		}
		return result;
	}
}

/*
2D Shape parallel to XZ plane
Can be used to plot 2D points with callback func
Return a 2D planer shape
*/
class Shape extends Origin{
	constructor(center, translation, rotation, callback){
		super(translation, rotation);
		this.center = center;
		this.callback = callback;
		this.pointList = [];
		this.segList = [];
		this.vertices2D = [];
		this.segList2D = [];
	}
	get centerPoint(){
		return new SinglePoint(this.center, this.translation, this.rotation);
	}
	//overloads x, y or x with function(X)
	pointTo(x, y){
		let t = 0;
		if (arguments.length === 1){
			if (this.callback && typeof this.callback === "function"){
				t = this.callback(x);
			}
			else {
				return "not valid function"
			}
		}
		else {
			t = y;
		}
		let newvector = [this.center[0]+x, this.center[1], this.center[2]+t]
		let apoint = new SinglePoint(newvector, this.translation, this.rotation)
		let prepoint = [];
		if (this.vertices2D.length !== 0){
			let x0 = this.vertices2D[this.vertices2D.length -1][0];
			let y0 = this.vertices2D[this.vertices2D.length -1][1];
			prepoint = [x0,y0];
			this.segList2D.push({a:prepoint, b:[x,t]})
			if (x=== this.vertices2D[0][0] && t === this.vertices2D[0][1]){
				return "shape closed"
			}
			else{
				this.pointList.push(apoint);
				this.vertices2D.push([x,t])
				return "point added"
			}
		}
		else {
			this.pointList.push(apoint);
			this.vertices2D.push([x,t])
			return "point added"
		}
	}
	//take out collinear points
	get excludeCollinear(){
		let list = this.segList2D;
		let newlist = [];
		for (let i=0; i<list.length; i++){
			let thisseg = list[i];
			let x0 = thisseg.a[0];
			let y0 = thisseg.a[1];
			let x1 = thisseg.b[0];
			let y1 = thisseg.b[1];
			let nextseg
			if (i === list.length-1){
				nextseg = list[0];
			}
			else{
				nextseg = list[i+1];
			}
			let a0 = nextseg.a[0];
			let b0 = nextseg.a[1];
			let a1 = nextseg.b[0];
			let b1 = nextseg.b[1];
			let xy = Math.round((y1-y0)/(x1-x0)*1000000)/1000000;
			let ab = Math.round((b1-b0)/(a1-a0)*1000000)/1000000;
			if (xy !== ab){
				newlist.push(thisseg);
			}
			else{
				if (i===list.length-1){
					newlist.shift();
				}
				newlist.push({a:[x0,y0],b:[a1,b1]})
			}
		}
		return newlist;
	}
	//clipp out two segments given the segment list
	clipEar(listofseg){
		let copypoint = listofseg.map(function(e){
			return e.a
		})
		for (let i=0; i<listofseg.length;i++){
			let sega = listofseg[i];
			let segb;
			let copy = copypoint.map((x=>x));
			if (i=== listofseg.length-1){
				segb = listofseg[0]
				copy.splice(i,1)
				copy.splice(0,1)
				copy.splice(0,1)
			}
			else {
				segb = listofseg[i+1]
				copy.splice(i,1)
				copy.splice(i,1)
				copy.splice(i,1)
			}
			let notfound = true;
			for (let j=0; j<copy.length;j++){
				if (this.checkTriangle(sega, segb, copy[j])){
					notfound = false
					break;
				}
			}
			if (notfound){
				this.segList.push([sega, segb])
				if (i=== listofseg.length-1){
					listofseg.pop();
					listofseg.shift();
					return listofseg.splice(0,0,{a:sega.a,b:segb.b})
				}
				else {
					return listofseg.splice(i,2,{a:sega.a,b:segb.b})
				}
			}
		}
	}
	get arrayBuffer(){
		let array = this.excludeCollinear
		if (array.length >= 4){
			do{
				this.clipEar(array)
			}
			while( array.length > 3 )
			this.segList.push([array[0],array[1]]);
			let result = new Float32Array(this.segList.length*9)
			let centerx = this.center[0];
			let centery = this.center[1];
			let centerz = this.center[2];
			let tran = this.translation;
			let rota = this.rotation;
			this.segList.forEach(function(e,i){
				let sega = e[0];
				let vectora = [sega.a[0] + centerx,centery,sega.a[1]+centerz]
				let vectorb = [sega.b[0] + centerx,centery,sega.b[1]+centerz]
				let segb = e[1];
				let vectorc = [segb.b[0] + centerx,centery,segb.b[1]+centerz]
				let tria = new Tria(vectora, vectorb, vectorc, tran, rota)
				result.set(tria.arrayBuffer,i*9)
			})
			this.segList = [];
			return result
		}
		else {
			throw array.length.toString() + " points are not enough"
		}
	}

	//check whether this point is inside of the triangle
	checkTriangle(x,y,z){
		let segmentA = x;
		let segmentB = y;
		let pointX = z;
		let pointA = segmentA.a
		let pointB = segmentA.b
		let pointC = segmentB.b
		let Ax = pointA[0];
		let Ay = pointA[1];
		let Bx = pointB[0];
		let By = pointB[1];
		let Cx = pointC[0];
		let Cy = pointC[1];
		let Xx = pointX[0];
		let Xy = pointX[1];
		let c = Math.sqrt((Ax-Bx)*(Ax-Bx)+(Ay-By)*(Ay-By))
		let a = Math.sqrt((Bx-Cx)*(Bx-Cx)+(By-Cy)*(By-Cy))
		let b = Math.sqrt((Ax-Cx)*(Ax-Cx)+(Ay-Cy)*(Ay-Cy))
		let xa = Math.sqrt((Xx-Ax)*(Xx-Ax)+(Xy-Ay)*(Xy-Ay))
		let xb = Math.sqrt((Xx-Bx)*(Xx-Bx)+(Xy-By)*(Xy-By))
		let xc = Math.sqrt((Xx-Cx)*(Xx-Cx)+(Xy-Cy)*(Xy-Cy))
		let gamma1 = Math.acos((xa*xa + xb*xb -c*c)/(2*xa*xb))
		let gamma2 = Math.acos((xb*xb + xc*xc -a*a)/(2*xb*xc))
		let gamma3 = Math.acos((xa*xa + xc*xc -b*b)/(2*xa*xc))
		let result = gamma1+gamma2+gamma3
		//on edge or point will return true
		if (Math.round((result - Math.PI*2)*1000)){
			return false;
		}
		else{
			return true;
		}
	}
}

/*
Spline is consisted of some segments
spline can take function, spline can be used by itself or used by extrusion
*/
class Spline extends Origin {
	constructor (center, t1, t2, callbackX, callbackY, callbackZ, segment, translation, rotation){
		super(translation, rotation);
		this.center = center;
		this.t1 = t1;
		this.t2 = t2;
		this.xfunction = callbackX;
		this.yfunction = callbackY;
		this.zfunction = callbackZ;
		this.segment = segment;
	}
	get segmentList(){
		let result = [];
		for (let i = 0; i< this.segment; i++){
			let ele = this.t1 + (this.t2- this.t1)/this.segment*i;
			let elex = this.center[0] + this.xfunction(ele);
			let eley = this.center[1] + this.yfunction(ele);
			let elez = this.center[2] + this.zfunction(ele);
			let vecta = [elex, eley, elez];
			let ele1 = this.t1 + (this.t2- this.t1)/this.segment*(i+1);
			let elex1 = this.center[0] + this.xfunction(ele1); 
			let eley1 = this.center[1] + this.yfunction(ele1);
			let elez1 = this.center[2] + this.zfunction(ele1);
			let vectb = [elex1, eley1, elez1]
			let aseg = new Segment(vecta, vectb, this.translation, this.rotation);
			result.push(aseg);
		}
		return result;
	}
	get arrayBuffer(){
		let result = new Float32Array(this.segmentList.length * 6);
		this.segmentList.forEach(function (e,i){
			let ele = e.arrayBuffer;
			result.set(ele,i*6);
		})
		return result;
	}
}

/*
Extrusion  will take a shape and extrude it along a spline
Shape's center point will become the control point of the extrusion
Each segment's starting point becomes the center point of the shape
*/
class Extrusion extends Spline {
	constructor (shape, spline){
		super(spline.center, spline.t1, spline.t2, spline.xfunction, spline.yfunction, 
		spline.zfunction, spline.segment, spline.translation, spline.rotation);
		this.shape = shape
	}
	get shapeList(){
		let result = [];
		for (let i=0; i<this.segmentList.length; i++){
			let seg = this.segmentList[i];
			let v1 = seg.vectora.map((x=>x));
			let v2 = seg.vectorb.map((x=>x))
			let v3 = [v1[0], 1, v1[2]];
			let atria = new Tria(v1, v2, v3, [0, 0, 0], [0, 1, 0, 0, "A"])
			let perpend = atria.perpendicularVector
			let xz = Math.sqrt((v2[0]-v1[0])*(v2[0]-v1[0])+(v2[2]-v1[2])*(v2[2]-v1[2]))
			let angle = Math.atan(xz/(v2[1]-v1[1])) * 180 /Math.PI
			let rotation = [angle, perpend[0], perpend[1], perpend[2],"A"]
			let newshape = new Shape(this.center, v1 , rotation)
			newshape.segList2D = this.shape.segList2D.map((x=>x))
			result.push(newshape)
		}
		return result;
	}
	get arrayBuffer(){
		let a = this.shapeList.length
		let b = this.shapeList[1].arrayBuffer.length
		let result = new Float32Array(a*b);
		this.shapeList.forEach(function(e,i){
			let ele = e.arrayBuffer
			result.set(ele, i*b)
		})
		return result
	}
	
	
}
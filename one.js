/*
One is designed to programatically create 3D points.
The end result is one dimentional array buffer.
*/

/*
Origin is the base class 
*/
class Origin {
	constructor (x, y, z){
		this.originX = x;
		this.originY = y;	
		this.originZ = z;	
	};
};

/*
Single Point is single vertex
array buffer as x, y, z
*/

class SinglePoint extends Origin {
	constructor (x, y, z, originX, originY, originZ) {
		super(originX, originY, originZ)
		this.X = x;
		this.Y = y;
		this.Z = z;
	}
	get arrayBuffer(){
		return [this.X, this.Y, this.Z];	
	} 
	get translationVector() {
		return [this.X - this.originX, this.Y- this.originY, this.Z- this.originZ];	
	} 
	
	get distance() {
		var xr = (this.X - this.originX) * (this.X - this.originX);
		var yr = (this.Y - this.originY) * (this.Y - this.originY);
		var zr = (this.Z - this.originZ) * (this.Z - this.originZ);	
		return Math.sqrt(xr+yr+zr);
	}
}


/*
Segment is consisted of two points. If two points are same return undefined
array buffer return A, B
*/

class Segment extends Origin {
	constructor (pointa, pointb, originX, originY, originZ){
		if (pointa[0] === pointb[0] && pointa[1] === pointb[1] && pointa[2] === pointb[2]){
			return;
		}
		else {
			super(originX, originY, originZ);
			this.pointA = new SinglePoint(pointa[0], pointa[1], pointa[2], originX, originY, originZ);
			this.pointB = new SinglePoint(pointb[0], pointb[1], pointb[2], originX, originY, originZ);
		}
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
		return this.pointA.translationVector; 
	}
	get vectorB(){
		return this.pointB.translationVector; 
	}
	collinear (x, y, z){
		let x0 = this.pointA.X;
		let y0 = this.pointA.Y; 
		let z0 = this.pointA.Z;
		let a = this.segmentVector[0];
		let b = this.segmentVector[1];
		let c = this.segmentVector[2];
		let t1, t2, t3
		//infinity does exist
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
	constructor (pointa, pointb, pointc, originx, originy, originz){
		if ((pointa[0] === pointb[0] && pointa[1] === pointb[1] && pointa[2] === pointb[2]) ||
			(pointa[0] === pointc[0] && pointa[1] === pointc[1] && pointa[2] === pointc[2]) ||
			(pointc[0] === pointb[0] && pointc[1] === pointb[1] && pointc[2] === pointb[2])){
			return;
		}
		let tempseg = new Segment(pointa, pointb, originx, originy, originz)
		if (tempseg.collinear(pointc[0], pointc[1], pointc[2])){
			return;
		}
		else {
			super(originx, originy, originz);
			this.pointA = new SinglePoint(pointa[0], pointa[1], pointa[2], originx, originy, originz); 
			this.pointB = new SinglePoint(pointb[0], pointb[1], pointb[2], originx, originy, originz);
			this.pointC = new SinglePoint(pointc[0], pointc[1], pointc[2], originx, originy, originz);	
		}
	}
	get segmentBC(){
		let bc = new Segment([this.pointB.X, this.pointB.Y, this.pointB.Z], [this.pointC.X, this.pointC.Y, this.pointC.Z], this.originX, this.originY, this.originZ);
		return bc;
	}
	get segmentAB(){
		let ab = new Segment([this.pointA.X, this.pointA.Y, this.pointA.Z], [this.pointB.X, this.pointB.Y, this.pointB.Z], this.originX, this.originY, this.originZ);
		return ab;
	}
	get segmentCA(){
		let ca = new Segment([this.pointC.X, this.pointC.Y, this.pointC.Z], [this.pointA.X, this.pointA.Y, this.pointA.Z], this.originX, this.originY, this.originZ);
		return ca;
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
	constructor (pointa, pointb, pointc, pointd, originx, originy, originz){
		if ((pointa[0] === pointb[0] && pointa[1] === pointb[1] && pointa[2] === pointb[2]) ||
			(pointa[0] === pointc[0] && pointa[1] === pointc[1] && pointa[2] === pointc[2]) ||
			(pointc[0] === pointb[0] && pointc[1] === pointb[1] && pointc[2] === pointb[2]) ||
			(pointc[0] === pointd[0] && pointc[1] === pointd[1] && pointc[2] === pointd[2]) ||
			(pointa[0] === pointd[0] && pointa[1] === pointd[1] && pointa[2] === pointd[2]) ||
			(pointb[0] === pointd[0] && pointb[1] === pointd[1] && pointb[2] === pointd[2])
			){
			return;
		}
		let temptri = new Tria(pointa, pointb, pointc, originx, originy, originz)
		if (temptri.coplanar(pointd[0], pointd[1], pointd[2])){
			super(pointa, pointb, pointc, originx, originy, originz);
			this.pointD = new SinglePoint(pointd[0], pointd[1], pointd[2], originx, originy, originz);
		}
		else {
			return;
		}
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z, this.pointC.X, this.pointC.Y, this.pointC.Z, 
		this.pointD.X, this.pointD.Y, this.pointD.Z]
	}
	get segmentCD(){
		let cd = new Segment(this.pointC, this.pointD, this.originX, this.originY, this.originZ)
		return cd;
	}
	get segmentDA(){
		let da = new Segment(this.pointD, this.pointA, this.originX, this.originY, this.originZ)
		return da;
	}
	
}


/*
One is designed to programatically create 3D points.
The end result is one dimentional array buffer.
*/

/*
Origin is the base class any 3D object
x, y, z translates, i, j, k rotates
*/
class Origin {
	constructor (x, y, z, i, j, k){
		this.originX = x;
		this.originY = y;
		this.originZ = z;
        this.originI = i;
        this.originJ = j;
        this.originK = k;
	};
};

/*
Single Point is single vertex
array buffer as x, y, z
*/

class SinglePoint extends Origin {
	constructor (x, y, z, originX, originY, originZ, i, j, k) {
		super(originX, originY, originZ, i, j, k);
        //this has to be rotated
        this.x = x;
		this.y = y;
		this.z = z;
	}
	//this has to be rotated
    get X(){
        return
    }
    get Y(){
        return
    }
    get Z(){
        return
    }
    
    get arrayBuffer(){
		return [this.X, this.Y, this.Z];
	}
    //this has to be rotated
	get translationVector() {
		return [this.X - this.originX, this.Y- this.originY, this.Z- this.originZ];
	}
	get originDistance() {
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
	constructor (vectora, vectorb, originX, originY, originZ, i, j, k){
		if (vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]){
			return;
		}
		else {
			super(originX, originY, originZ, i, j, k);
			this.pointA = new SinglePoint(vectora[0], vectora[1], vectora[2], originX, originY, originZ, i, j, k);
			this.pointB = new SinglePoint(vectorb[0], vectorb[1], vectorb[2], originX, originY, originZ, i, j, k);
		}
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z];
	}
	//this needs to be rotated
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
	constructor (vectora, vectorb, vectorc, originx, originy, originz, i, j, k){
		if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2])){
			return;
		}
		let tempseg = new Segment(vectora, vectorb, 0, 0, 0, 0, 0, 0)
		if (tempseg.collinear(vectorc[0], vectorc[1], vectorc[2])){
			return;
		}
		else {
			super(originx, originy, originz, i, j, k);
			this.pointA = new SinglePoint(vectora[0], vectora[1], vectora[2], originx, originy, originz, i, j, k);
			this.pointB = new SinglePoint(vectorb[0], vectorb[1], vectorb[2], originx, originy, originz, i, j, k);
			this.pointC = new SinglePoint(vectorc[0], vectorc[1], vectorc[2], originx, originy, originz, i, j, k);
		}
	}
	get segmentBC(){
		let bc = new Segment(this.pointB.translationVector, this.pointC.translationVector, this.originX, this.originY, this.originZ, this.originI, this.originJ, this.originK);
		return bc;
	}
	get segmentAB(){
		let ab = new Segment(this.pointA.translationVector, this.pointB.translationVector, this.originX, this.originY, this.originZ, this.originI, this.originJ, this.originK);
		return ab;
	}
	get segmentCA(){
		let ca = new Segment(this.pointC.translationVector, this.pointC.translationVector, this.originX, this.originY, this.originZ, this.originI, this.originJ, this.originK);
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
	constructor (vectora, vectorb, vectorc, vectord, originx, originy, originz, i, j, k){
		if ((vectora[0] === vectorb[0] && vectora[1] === vectorb[1] && vectora[2] === vectorb[2]) ||
			(vectora[0] === vectorc[0] && vectora[1] === vectorc[1] && vectora[2] === vectorc[2]) ||
			(vectorc[0] === vectorb[0] && vectorc[1] === vectorb[1] && vectorc[2] === vectorb[2]) ||
			(vectorc[0] === vectord[0] && vectorc[1] === vectord[1] && vectorc[2] === vectord[2]) ||
			(vectora[0] === vectord[0] && vectora[1] === vectord[1] && vectora[2] === vectord[2]) ||
			(vectorb[0] === vectord[0] && vectorb[1] === vectord[1] && vectorb[2] === vectord[2])
			){
			return;
		}
		let temptri = new Tria(vectora, vectorb, vectorc, originx, originy, originz, i, j, k)
		if (temptri.coplanar(vectord[0]+originx, vectord[1]+originy, vectord[2]+originz)){
			super(vectora, vectorb, vectorc, originx, originy, originz);
			this.pointD = new SinglePoint(vectord[0], vectord[1], vectord[2], originx, originy, originz);
		}
		else {
			return;
		}
	}
	get arrayBuffer(){
		return [this.pointA.X, this.pointA.Y, this.pointA.Z, this.pointB.X, this.pointB.Y, this.pointB.Z, this.pointC.X, this.pointC.Y, this.pointC.Z, 
		this.pointD.X, this.pointD.Y, this.pointD.Z];
	}
	get segmentCD(){
		let cd = new Segment(this.pointC.translationVector, this.pointD.translationVector, this.originX, this.originY, this.originZ)
		return cd;
	}
	get segmentDA(){
		let da = new Segment(this.pointD.translationVector, this.pointA.translationVector, this.originX, this.originY, this.originZ)
		return da;
	}
}

/*
Orthogonal Box is consisted of two Quads
Origin stays in center 
*/

class Box extends Origin {
	constructor (length, width, height, originx, originy, originz, i, j, k){
		super(originx, originy, originz, i, j, k);
		this.length = length;
		this.width = width;
		this.height = height;
	}
	get topQuad(){
		let A = [- this.width/2, this.height/2, - this.length/2];
		let B = [this.width/2, this.height/2, - this.length/2];
		let C = [this.width/2, this.height/2, this.length/2];
		let D = [- this.width/2, this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.originX, this.originY, this.originZ, this.I, this.J, this.K);
		return result;
	}
	get bottomQuad(){
		let A = [- this.width/2, - this.height/2, - this.length/2];
		let B = [this.width/2, - this.height/2, - this.length/2];
		let C = [this.width/2, - this.height/2, this.length/2];
		let D = [- this.width/2, - this.height/2, this.length/2];
		let result = new Quad(A,B,C,D,this.originX, this.originY, this.originZ, this.I, this.J, this.K);
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
	constructor (edge, radius, centervector, normalvector, originx, originy, originz){
		super(originx, originy, originz);
		this.edge = edge;
		this.radius = radius;
        this.center = centervector;
        this.normal = normalvector;
	}
	getAxisIntersect(axis) {
		let a = this.originX;
		let b = this.originY;
		let c = this.originZ;
		switch (axis){
			case "X":
			return (c*c + b*b)/a + a;
			break;
			case "Y":
			return (c*c + a*a)/b + b;
			break;
			case "Z":
			return (b*b + a*a)/c + c;	
		}
	}
	get triaList() {
		let result = [];
		let xinter = this.getAxisIntersect("X");
		let yinter = this.getAxisIntersect("Y");
		let zinter = this.getAxisIntersect("Z");
		//only intersect with Y axis
        if (xinter === Infinity && zinter === Infinity && yinter !== Infinity) {
            for (let i=0; i<this.edge; i++){
                let ele = new Tria([0, this.originY,0], [this.radius, this.originY, 0])     
                
            }
            
        }
        //only intersect with X axis
        else if(yinter === Infinity && zinter === Infinity && xinter !== Infinity){
            
        }
        //only intersect with Z axis
        else if(xinter === Infinity && yinter === Infinity && zinter !== Infinity){
            
        }
        //intersect with XY axis
        else if(xinter !== Infinity && yinter !== Infinity && zinter === Infinity){
            
        }
        //intersect with YZ axis
        else if(yinter !== Infinity && zinter !== Infinity && xinter === Infinity){
            
        }
        //intersect with XZ axis
        else if(xinter !== Infinity && zinter !== Infinity && yinter === Infinity){
            
        }
        //intersect with XYZ axis
        else {
            
        }
	}
	
}

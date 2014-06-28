var D=2; // There is always an extra 'height' dimension here
var m=2; // Minimum height difference for avalanche
var N=4; // Usual molecule linear dimension

var moleculeHash = [];
var transferMatrix = [];
var transferCount = 0;
var transferSubtotal = [];

function clone(obj) {
    var c = {};
    for(var i in obj) {
        if(typeof(obj[i])=="object" && obj[i] != null)
            c[i] = clone(obj[i]);
        else
            c[i] = obj[i];
    }
    return c;
}

function Molecule(lengtharray)
{
	this.changed = 1;
	this.n = lengtharray;
	this.sandArray = [];	
	this.donor = [];
	this.acceptor = [];
	
	this.totalSize = 1;
	
	for (var i=0;i<D;i++)
		this.totalSize *= lengtharray[i];
	
	for (var i=0;i<this.totalSize;i++)
	{
		this.sandArray[i]=0;	
		this.donor[i] = -1;
		this.acceptor[i] = -1;
	}
}

Molecule.prototype.toHash = function()
{
	var hashStr = "";
	
	for (var i=0;i<this.totalSize;i++)
		hashStr = hashStr + this.sandArray[i] + ",";
		
	return hashStr;
}

Molecule.prototype.coordToIndex = function(coord)
{
	var idx = 0;
	var modulator = 1;
	
	for (var i=0;i<D;i++)
	{
		idx += coord[i] * modulator;
		modulator *= this.n[i];
	}
	
	return idx;
}

Molecule.prototype.indexToCoord = function(idx)
{
	var coord = [];
	var divisor = 1;
	
	for (var i=0;i<D;i++)
	{		
		coord[i] = (Math.floor(idx/divisor))%this.n[i];
		divisor *= this.n[i];
	}
	
	return coord;
}

Molecule.prototype.getCell = function(coord)
{
	var idx = this.coordToIndex(coord);
	return this.sandArray[idx];
}

Molecule.prototype.setCell = function(coord, val)
{
	var idx = this.coordToIndex(coord);
	this.sandArray[idx] = val;
}

Molecule.prototype.join = function(othermolecule, edgeidx)
{
	var newarray = clone(this.n);
	
	newarray[edgeidx] += othermolecule.n[edgeidx];
	
	var newmolecule = new Molecule(newarray);
	
	for (var i=0;i<newmolecule.totalSize;i++)
	{
		var x=newmolecule.indexToCoord(i);
		
		if (x[edgeidx] < this.n[edgeidx])
		{
			newmolecule.setCell(x, this.getCell(x));
		}
		else
		{
			var y = clone(x);
			
			y[edgeidx] -= this.n[edgeidx];
			newmolecule.setCell(x, othermolecule.getCell(y));
		}
	}
	
	return newmolecule;
}

Molecule.prototype.split = function(edgeidx)
{
	var newarray = clone(this.n);
	newarray[edgeidx] /= 2;
	
	var newmolecule1 = new Molecule(newarray);
	var newmolecule2 = new Molecule(newarray);
	
	for (var i=0;i<newmolecule1.totalSize;i++)
	{
		var x = newmolecule1.indexToCoord(i);
		var y = clone(x);
		
		y[edgeidx] += newmolecule1.n[edgeidx];
		
		newmolecule1.setCell(x, this.getCell(y));
		newmolecule2.setCell(x, this.getCell(x));
	}
	
	return [newmolecule1, newmolecule2];
}

Molecule.prototype.getAdjacentSites = function(idx)
{
	var coord = this.indexToCoord(idx);
	var sitelist = [];
	
	for (var i=0;i<D;i++)
	{		
		coord[i]++; if (coord[i] < this.n[i]) sitelist.push(this.coordToIndex(coord)); coord[i]--;
		coord[i]--; if (coord[i] >= 0) sitelist.push(this.coordToIndex(coord)); coord[i]++;
	}
	
	return sitelist;
}

Molecule.prototype.getDonors = function()
{
	for (var i=0;i<this.totalSize;i++)
	{
		var sitelist = this.getAdjacentSites(i);
		
		var tie = 0;
		var bestheight = 0, bestidx = -1;
		var height0 = this.sandArray[i];
		
		for (var j=0;j<sitelist.length;j++)
		{
			var delta = this.sandArray[sitelist[j]] - height0;
			
			if (delta >= m)			
			{
				if (delta == bestheight)
				{
					tie = 1;
				}
				else if (delta > bestheight)
				{
					tie = 0;
					bestheight = delta;
					bestidx = sitelist[j];
				}
			}
		}
		
		if ((!tie)&&(bestidx>=0))
			this.donor[i] = bestidx;
		else this.donor[i] = -1;
	}
}

Molecule.prototype.getAcceptors = function()
{
	for (var i=0;i<this.totalSize;i++)
	{
		var sitelist = this.getAdjacentSites(i);
		var acceptorList = [];
		
		for (var j=0;j<sitelist.length;j++)
			if (this.donor[sitelist[j]]==i)
				acceptorList.push(sitelist[j]);
				
		var tie = 0;
		var bestheight = 0, bestidx = -1;
		var height0 = this.sandArray[i];
		
		for (var j=0;j<acceptorList.length;j++)
		{
			var delta = height0 - this.sandArray[acceptorList[j]];
			
			if (delta == bestheight)
			{
				tie = 1;
			}
			else if (delta > bestheight)
			{
				tie = 0;
				bestheight = delta;
				bestidx = acceptorList[j];
			}
		}
		
		if ((!tie)&&(bestidx>=0))
			this.acceptor[i] = bestidx;
		else this.acceptor[i] = -1;
	}
}

Molecule.prototype.iterateAvalanche = function()
{
	var changes = 0;
	
	this.getDonors();
	this.getAcceptors();
	
	for (var i=0;i<this.totalSize;i++)
	{
		if (this.acceptor[i]>=0)
		{
			changes++;
			
			this.sandArray[i]--;
			this.sandArray[this.acceptor[i]]++;
		}
	}
	
	return changes;
}

Molecule.prototype.getReactivityParams = function()
{
	var hbar=0, hmin=10000, hmax=0, h2bar=0, count=0;
	
	for (var i=0;i<this.totalSize;i++)
	{
		hbar += this.sandArray[i];
		h2bar += this.sandArray[i]*this.sandArray[i];
		if (this.sandArray[i] < hmin) hmin = this.sandArray[i];
		if (this.sandArray[i] > hmax) hmax = this.sandArray[i];
		count++;
	}
	
	h2bar /= count; hbar /= count;
	
	return { average: hbar, stdev: Math.sqrt(h2bar-hbar*hbar), min: hmin, max: hmax };
}

var boundaryMolecules = [];
var Bath = [];

function initBoundaries()
{
	var molSize = [];
	for (var i=0;i<D;i++)
	{
		molSize[i] = N;
	}
	
	boundaryMolecules.push(new Molecule(molSize));
	boundaryMolecules.push(new Molecule(molSize));
	
	for (var i=0;i<boundaryMolecules[0].totalSize;i++)
		boundaryMolecules[0].sandArray[i] = Math.floor(3*Math.random()); 

	for (var i=0;i<boundaryMolecules[1].totalSize;i++)
		boundaryMolecules[1].sandArray[i] = 8+Math.floor(3*Math.random()); // This one starts a bit higher up
}	

function iterateBath()
{
	var len0=Bath.length;
	var subBath = [];
	
	// Iterate each molecule in isolation 4 times
	for (var i=0;i<len0;i++)
	{
		if (Bath[i].changed)
		{
			do
			{
				Bath[i].changed = Bath[i].iterateAvalanche();
			} while (Bath[i].changed);
		}
	}
	
	// Pick some pairs, merge, iterate, and split
	if (len0>1)
	for (var i=0;i<len0/4;i++)
	{
		var m1 = Math.floor(Math.random()*Bath.length);
		var m2;
		
		do
		{
			m2 = Math.floor(Math.random()*Bath.length);
		} while (m2==m1);
		
		var dir = Math.floor(Math.random()*D);
		var joinedMolecule = Bath[m1].join(Bath[m2], dir);

		do
		{
			joinedMolecule.changed = joinedMolecule.iterateAvalanche();
		} while (joinedMolecule.changed);
		
		var products = joinedMolecule.split(dir);
		
		var rparam_old = [ Bath[m1].getReactivityParams(), Bath[m2].getReactivityParams() ];		
		var rparam_new = [ products[1].getReactivityParams(), products[0].getReactivityParams() ];
		
		idx1 = Math.floor(rparam_old[0].stdev*50/4);
		idx2 = Math.floor(rparam_new[0].stdev*50/4);
		
		if ((idx1>=0)&&(idx2>=0)&&(idx1<50)&&(idx2<50))
		{
			transferSubtotal[idx1]++;
			
			transferMatrix[idx1+idx2*50].y++;
		}
		transferCount++;

		idx1 = Math.floor(rparam_old[1].stdev*50/4);
		idx2 = Math.floor(rparam_new[1].stdev*50/4);
		
		if ((idx1>=0)&&(idx2>=0)&&(idx1<50)&&(idx2<50))
		{
			transferSubtotal[idx1]++;
			transferMatrix[idx1+idx2*50].y++;
		}
		transferCount++;
		
		Bath[m1]=products[1];
		Bath[m2]=products[0];
		Bath[m1].changed=1; Bath[m2].changed=1;
	}
	
	// Also do it with boundary molecules
	for (var i=0;i<len0/4;i++)
	{
		var joinedMolecule;
		var m1 = Math.floor(Math.random()*Bath.length);
		var m2 = Math.floor(Math.random()*2);
		
		if (m2 == 0)
		{			
			joinedMolecule = Bath[m1].join(boundaryMolecules[0],0);
			do
			{
				joinedMolecule.changed = joinedMolecule.iterateAvalanche();
			} while (joinedMolecule.changed);
			
			var products = joinedMolecule.split(0);
		
			var rparam_old = Bath[m1].getReactivityParams();		
			var rparam_new = products[1].getReactivityParams();
			
			idx1 = Math.floor(rparam_old.stdev*50/4);
			idx2 = Math.floor(rparam_new.stdev*50/4);
		
			if ((idx1>=0)&&(idx2>=0)&&(idx1<50)&&(idx2<50))
			{
				transferSubtotal[idx1]++;
				transferMatrix[idx1+idx2*50].y++;
			}
			transferCount++;
		
			Bath[m1] = products[1];
		}
		else
		{
			joinedMolecule = boundaryMolecules[1].join(Bath[m1],0);
			do
			{
				joinedMolecule.changed = joinedMolecule.iterateAvalanche();
			} while (joinedMolecule.changed);

			var products = joinedMolecule.split(0);
			var rparam_old = Bath[m1].getReactivityParams();		
			var rparam_new = products[0].getReactivityParams();
			
			idx1 = Math.floor(rparam_old.stdev*50/4);
			idx2 = Math.floor(rparam_new.stdev*50/4);
		
			if ((idx1>=0)&&(idx2>=0)&&(idx1<50)&&(idx2<50))
			{
				transferSubtotal[idx1]++;
				transferMatrix[idx1+idx2*50].y++;
			}
			transferCount++;

			Bath[m1] = products[0];
		}
		
		Bath[m1].changed=1;
	}
}

function initSandpile()
{
	var molSize = [];
	for (var i=0;i<D;i++)
	{
		molSize[i] = N;
	}
	
	for (var i=0;i<1000;i++)
	{
		var mol = new Molecule(molSize);
		
		for (var j=0;j<mol.totalSize;j++)
		{
			mol.sandArray[j] = Math.floor(Math.random()*10);
		}
		
		mol.changed = 1;
		
		Bath.push(mol);
	}
}

var height_dist = [];
var stdev_dist = [];
var dcount = 0;

function collectStatistics()
{	
	// Determine distribution of molecules
	for (var i=0;i<Bath.length;i++)
	{
		moleculeHash[Bath[i].toHash()]++;
		
		var rparams = Bath[i].getReactivityParams();
		
		var idx = Math.floor(rparams.average*50/10);
		if (idx<50) { height_dist[idx].y++; };
		
		idx = Math.floor(rparams.stdev*50/4);
		if (idx<50) { stdev_dist[idx].y++; };
		
		dcount++;
	}	
}

var graph = {};

var margin = {top: 20, right: 20, bottom: 40, left: 40},
		width = 250 - margin.left - margin.right,
		height = 250 - margin.top - margin.bottom;

function setupGraph()
{

	graph.x = d3.scale.linear()
				.range([0, width]);
	graph.x2 = d3.scale.linear()
				.range([0, width]);

	graph.y = d3.scale.linear()
				.range([height, 0]);
	graph.y2 = d3.scale.linear()
				.range([height, 0]);
				
	graph.xAxis = d3.svg.axis()
				    .scale(graph.x)
				    .orient("bottom");
	graph.xAxis2 = d3.svg.axis()
				    .scale(graph.x2)
				    .orient("bottom");
	graph.yAxis = d3.svg.axis()
				    .scale(graph.y)
				    .orient("left");
	graph.yAxis2 = d3.svg.axis()
				    .scale(graph.y2)
				    .orient("left");

	graph.x.domain([0,10]);
	graph.x2.domain([0,4]);
	graph.y.domain([0,2]);
	graph.y2.domain([0,3]);
	
	var g1 = d3.select("#svg1")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
		.append("g")
		.attr("class", "graph1")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
		
	var g2=d3.select("#svg2")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
		.append("g")
		.attr("class", "graph2")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
		
	
	g1
		.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + height + ")")
		.call(graph.xAxis)
		.append("text")
		.attr("transform", "translate(0, 30)")		
		.text("Average Height");

   g1
	  .append("g")
      .attr("class", "y axis")
      .call(graph.yAxis)
	  .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("PDF");      

	g2
		.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + height + ")")
		.call(graph.xAxis2)
		.append("text")
		.attr("transform", "translate(0, 30)")		
		.text("Height Stdev");
		
   g2
	  .append("g")
      .attr("class", "y axis")
      .call(graph.yAxis2)
	  .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("PDF");      

	zeroStatistics(); dcount = 1;
	
	var svg1 = d3.select(".graph1");
	var svg2 = d3.select(".graph2");
	svg1.selectAll(".bar")
		.data(height_dist)
		.enter().append("rect")
		.attr("class", "bar")
		.attr("x", function(d) { return graph.x(d.x); })
		.attr("width", 4)
		.attr("y", function(d) { return graph.y( (50/10) * d.y / dcount); })
		.attr("height", function(d) { return height - graph.y( (50/10) * d.y / dcount); });

	svg2.selectAll(".bar")
		.data(stdev_dist)
		.enter().append("rect")
		.attr("class", "bar")
		.attr("x", function(d) { return graph.x2(d.x); })
		.attr("width", 4)
		.attr("y", function(d) { return graph.y2( (50/4) * d.y / dcount); })
		.attr("height", function(d) { return height - graph.y2( (50/4) * d.y / dcount); });


	var g3 = d3.select("#svg3")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
		.append("g")
		.attr("class", "graph3")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	g3.selectAll(".bar")
		.data(transferMatrix)
		.enter()
		.append("rect")
		.attr("class","gsquare")
		.attr("x", function(d, i) { var x=i%50; return 200*x/50; })
		.attr("y", function(d, i) { var y=Math.floor(i/50)%50; return 200 - 200*y/50; })
		.attr("width", 5)
		.attr("height", 5)
		.attr("fill", function(d, i) { 
			var val = d.y;			
			if (transferSubtotal[i%50]>0) val/=transferSubtotal[i%50];
			if (val>1) val=1; if (val<0) val=0;
			return "rgb("+(255*val)+",0,"+(255*(1-val))+")"; 
		});
		
	g3.append("text")
		.text("Initial Stdev")
		.attr("transform", "translate(60,"+(height+35)+")");

	g3.append("text")
		.text("Resultant Stdev")
		.attr("transform", "rotate(-90), translate(-150, -10)");
}

function displayData()
{		
	var svg1 = d3.select(".graph1");
	var svg2 = d3.select(".graph2");
	var svg3 = d3.select(".graph3");
	
	svg1.selectAll(".bar")
		.data(height_dist)
		.attr("y", function(d) { return graph.y( (50/10) * d.y / dcount); })
		.attr("height", function(d) { return height - graph.y( (50/10) * d.y / dcount); });

	svg2.selectAll(".bar")
		.data(stdev_dist)
		.attr("y", function(d) { return graph.y2( (50/4) * d.y / dcount); })
		.attr("height", function(d) { return height - graph.y2( (50/4) * d.y / dcount); });	
		
	svg3.selectAll(".gsquare")
		.data(transferMatrix)
		.attr("fill", function(d, i) { 
			var val = d.y;
			var j = i%50; //Math.floor(i/50)%50;
			if (transferSubtotal[j]>0) val /= transferSubtotal[j];
			
			val = (Math.log(val+1e-4)/Math.log(10)+4)/4.0;
			
			if (val>1) val=1; if (val<0) val=0;
			
			var r = Math.floor(255*val);
			var b = Math.floor(255*(1-val));
			
			return "rgb("+r+",0,"+b+")"; 
		});
			
	molCount = 0;
	for (var key in moleculeHash)
	{
		molCount++;
	}
	document.getElementById("diversity").innerHTML = "" + molCount;
}

function zeroStatistics()
{
		
	for (var i=0;i<50;i++)
	{
		height_dist[i] = { x: i*10/50, y: 0 };
		stdev_dist[i] = { x: i*4/50, y: 0 };
	}
	
	dcount = 0;
}

function simulationStep()
{
	zeroStatistics();
//	for (var i=0;i<3;i++)
	{
		iterateBath();
		collectStatistics();
	}
	
	displayData();
	
	window.requestAnimationFrame(simulationStep);
}

window.onload = function()
{
	for (var i=0;i<50*50;i++)
		transferMatrix[i] = { y: 0 };
	for (var i=0;i<50;i++)
		transferSubtotal[i] = 0;
		
	setupGraph();
	initSandpile();
	initBoundaries();
	
	window.requestAnimationFrame(simulationStep);
}

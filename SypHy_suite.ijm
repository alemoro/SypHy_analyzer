// SypHy Analyzer

/*
This toolset implement different functionality for the analysis of SypHy.

Start developing 2018.01.25

Modify           
	2018.01.25 - first functionality: detect synapses and save them
	2018.02.05 - update ROI detection with DFF0
	2018.02.06 - add filter to ROI detection
*/

var majVer = 0;
var minVer = 15;
var about = "Developed by Alessandro Moro<br>"
			+ "<i>Department of Functional Genomics</i> (FGA)<br>"
			+ "<i>Centre of neuroscience and cognitive research</i> (CNCR)<br>"
			+ "<i>Vrij Universiteit</i> (VU) Amsterdam.<br>"
			+ "<i>email: a.moro@vu.nl</i><br><br><br>";

			
// Initialize the variables
var nh4Start = 200;
var sigma = 0.65;
var thrRadius = 7;
var minSize = 9;
var nStd = 5;

bFirst = call("ij.Prefs.get", "syphyOption.FirstRun", true);
if(bFirst){
	call("ij.Prefs.set", "syphyOption.FirstRun", false);
}else{
	nh4Start = call("ij.Prefs.get", "syphyOption.nh4Start", true);
	sigma = call("ij.Prefs.get", "syphyOption.sigma", true);
	thrRadius = call("ij.Prefs.get", "syphyOption.thrRadius", true);
	minSize = call("ij.Prefs.get", "syphyOption.minSize", true);
	nStd = call("ij.Prefs.get", "syphyOption.nStd", true);
}

// Start for now
Dialog.create("Set SypHy parameters");
	Dialog.addNumber("NH4 start", nh4Start);
	Dialog.addNumber("Gaussian sigma", sigma);
	Dialog.addNumber("Threshold radius", thrRadius);
	Dialog.addCheckbox("Process folder", false);
	Dialog.addNumber("Minimum synapse size", minSize);
	Dialog.addNumber("Number of STD", nStd);
	Dialog.show;
	nh4Start = Dialog.getNumber();
	sigma = Dialog.getNumber();
	thrRadius = Dialog.getNumber();
	bBatch = Dialog.getCheckbox();
	minSize = Dialog.getNumber();
	nStd = Dialog.getNumber();

call("ij.Prefs.set", "syphyOption.nh4Start", nh4Start);
call("ij.Prefs.set", "syphyOption.sigma", sigma);
call("ij.Prefs.set", "syphyOption.thrRadius", thrRadius);
call("ij.Prefs.set", "syphyOption.minSize", minSize);
call("ij.Prefs.set", "syphyOption.nStd", nStd);

if(bBatch){
	batchDir = getDirectory("Select SypHy folder");
	batchList = getFileList(batchDir);
	nBatch = batchList.length;
}else{
	nBatch = 1;
}

setBatchMode(true);
for(b=0; b<nBatch; b++){
	if(bBatch){
		if(endsWith(batchList[b], ".tif") || endsWith(batchList[b], ".stk")){
			open(batchDir + batchList[b]);
		}
	}else {
		setBatchMode("hide");
	}
	orTitle = getTitle();
	// First apply a DF/F0
	run("Z Project...", "stop=20 projection=[Average Intensity]");
	rename("F0");
	imageCalculator("Subtract create 32-bit stack", orTitle,"F0");
	rename("DF");
	imageCalculator("Divide create 32-bit stack", "DF","F0");
	rename("DFF0");
	selectWindow("F0");
	close();
	selectWindow("DF");
	close();
	selectWindow("DFF0");
	run("Gaussian Blur 3D...", "x=" + sigma +" y=" + sigma +" z=" + sigma);
	run("Z Project...", "start=" + nh4Start + " projection=[Max Intensity]");
	//run("Unsharp Mask...", "radius="+thrRadius+" mask=0.40");
	run("Duplicate...", " ");
	run("Gaussian Blur...", "sigma="+minSize);
	imageCalculator("Subtract", "MAX_DFF0","MAX_DFF0-1");
	selectWindow("MAX_DFF0-1"); close();
	getStatistics(dff0A, dff0Mean, dff0Min, dff0Max, dff0Std, dff0Hist);
	setMinAndMax(dff0Mean, dff0Max);
	run("Morphological Filters", "operation=[White Top Hat] element=Disk radius="+thrRadius-1);
	/*
	run("Duplicate...", "duplicate range="+ nh4Start-1);
	rename("tempDiff");
	run("Gaussian Blur 3D...", "x=" + sigma +" y=" + sigma +" z=" + sigma);
	setPasteMode("Subtract");
	run("Set Slice...", "slice="+nSlices);
	run("Select All");
	for(i=nSlices; i>1; i--) {
		setSlice(i-1);
		run("Copy");
	   	setSlice(i);
		run("Paste");
	}
	run("Z Project...", "start=2 projection=[Max Intensity]");
	*/
	rename("Mask");
	run("8-bit");
	run("Auto Local Threshold", "method=Bernsen radius=" + thrRadius + " parameter_1=0 parameter_2=0 white");
	//selectWindow("tempDiff");
	//close();
	selectWindow("Mask");
	run("Analyze Particles...", "size="+minSize+"-200 add");
	selectWindow("Mask");
	close();
	selectWindow("DFF0");
	close();
	selectWindow(orTitle);
	selectWindow("MAX_DFF0");
	close();
	r = 0;
	nRoi = roiManager("Count");
	while(r<nRoi){
		c=0;
		roiManager("Select", r);
		Zprof = newArray(nSlices);
		for(s=1;s<=nSlices;s++){
			setSlice(s);
			getStatistics(zArea, Zmean, Zmin, Zmax, Zstd, Zhist);
			Zprof[s-1] = Zmean;
		}
		Zprof = Array.slice(Zprof,0,200);
		// add a static std
		zBase = Array.slice(Zprof,0,20);
		Array.getStatistics(zBase, statMin, statMax, statMean, statStd);
		stMean = newArray(Zprof.length);
		stStd = newArray(Zprof.length);
		Array.fill(stMean, statMean);
		Array.fill(stStd, statMean+nStd*statStd);
		// check if it goes above
		elCount = newArray(stStd.length);
		for(e=0;e<stMean.length;e++){
			if(stStd[e]<Zprof[e]){
				cc++;
			} else {
				cc = 0;			
			}
			elCount[e] = cc;
		}
		Array.getStatistics(elCount, elMin, elMax, elMean, elStd);
		if(elMax<5){
			roiManager("Delete");
			nRoi = roiManager("Count");
		} else {
			r++;
		}
	}
	newImage("Untitled", "8-bit black", 512, 512, 1);
	roiManager("deselect");
	roiManager("Fill");
	run("Dilate");
	run("Watershed");
	roiManager("reset");	
	selectWindow(orTitle);
	run("Remove Overlay");
	selectWindow("Untitled");
	run("Analyze Particles...", "size="+minSize+"-200 add");
	selectWindow("Untitled");
	close();
	
	if(bBatch){
		wd = getDirectory("Image");
		roiManager("Save", wd + "RoiSet_" + substring(orTitle,0,lengthOf(orTitle)-4) + ".zip");
		roiManager("Reset");
		while(nImages>0){
			close();
		}
	}else {
		setBatchMode("exit and display");
	}
}
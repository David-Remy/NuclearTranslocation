/* This macro is designed to measure protein intensity in nucleus and in cytoplasm and calculates the ratio
 *  Works on ImageJ 1.53c
 * macro_NuclearTranslocation.ijm
 *  Crated by DR (Nov 2021)
 */

//Clear all 
setOption("ExpandableArrays", true);
roiManager("reset");
run("Clear Results");
run("Close All");
print("\\Clear");
run("Set Measurements...", "area mean standard modal min centroid integrated area_fraction redirect=None decimal=3");

/////////////////////////////////////////////////////////////
////// begining of parameters customizable by the user //////
/////////////////////////////////////////////////////////////
// Radius for gaussian blur filter on nucleis, in pixels
sigma = 2;
// minimum area for nuclei, in pixels
minSizeNuclei = 800;
/////////////////////////////////////////////////////////////
//////// end of parameters customizable by the user /////////
/////////////////////////////////////////////////////////////

// Select input directory and create a "Translocation_Analysis" output directory in which everything will be saved 
dir_input=getDirectory("Select input folder");
condition=File.getName(dir_input);

if( !File.exists(dir_input+"Translocation_Analysis") ) {
	File.makeDirectory(dir_input+"Translocation_Analysis");
}

dir_output=dir_input+"Translocation_Analysis"+File.separator; 

//Asks the user to enter parameters: the extension of the masterfile, and designate the channels of interest for the nucleus and for TFEB
projection=false;
file_extensions=newArray(".nd2", ".nd", ".tif");
channels_array=newArray(1, 2, 3, 4);

Dialog.createNonBlocking("Channels");
Dialog.addChoice("Extension of files", file_extensions);
Dialog.addMessage("Channels ?");
Dialog.addMessage("1: DAPI ; 2:GFP ; 3: Cy3; 4: Cy5"); 
Dialog.addChoice("Nucleus channel", channels_array, "1.0");
Dialog.addChoice("Cell contour channel", channels_array, "3.0");
Dialog.addChoice("Translocated protein channel", channels_array, "2.0");
Dialog.addCheckbox("Images are Z-Stacks", true);

Dialog.show();

extension=Dialog.getChoice();
Nuc_chan=Dialog.getChoice();
contour_chan=Dialog.getChoice();
TrPr_chan=Dialog.getChoice();
use_project=Dialog.getCheckbox();
 
//Gets the list of all the files in the input directory
Filelist=getFileList(dir_input);
Array.sort(Filelist);
run("Set Measurements...", "area mean standard modal min centroid integrated area_fraction redirect=None decimal=3");

//Initializes arrays to compile all names of images, intensity of nuclei and of cytoplasm.
Name=newArray();
CellArea=newArray();
NucleiMean=newArray();
TotalMean=newArray();

// Open each acquisition and measure intensity for nuclei and cytoplasm 

for (i_files=0; i_files<lengthOf(Filelist); i_files++) {
	if (endsWith(Filelist[i_files], extension)) { // Opens .nd or .nd2 or .tif 
		name=File.getName(dir_input+Filelist[i_files]);
		subname=substring(name, 0, lastIndexOf(name, "."));	
		run("Bio-Formats", "open=["+dir_input+Filelist[i_files]+"] color_mode=Composite open_files view=Hyperstack stack_order=XYCZT");
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); 
		getDimensions(width, height, channels, slices, frames);
		
		if (slices>1) {
		
			Stack.setChannel(TrPr_chan);
			run("Enhance Contrast", "saturated=0.35");
			Dialog.createNonBlocking("Plans choice");
			Dialog.addMessage("Choose the plan to analyze, between 1 and "+slices+" (put the same value if you want only one image)."); 
			Dialog.addMessage("If several plans, a projection will be performed.");
			Dialog.addNumber("First plan (above 1)", 1);
			Dialog.addNumber("Last Plan (below "+slices+")", slices);
			Dialog.show();
		
			first_plan = Dialog.getNumber();
			last_plan = Dialog.getNumber();
			
			print(first_plan);
			print(last_plan);
			
			if( first_plan < 1 || first_plan > slices || last_plan < 1 || last_plan > slices) {
				exit("You choose an uncorrect number of plan");
			}
			
			
		
			selectWindow(name);
			run("Z Project...", "start="+first_plan+" stop="+last_plan+" projection=[Max Intensity]");
			close(name);
			selectWindow("MAX_"+name);
			rename(name);
		}
		
		selectWindow(name);
		run("Duplicate...", "title=Translocated_Prot duplicate channels="+TrPr_chan);
		selectWindow(name);
		run("Duplicate...", "title=nucleus duplicate channels="+Nuc_chan);
		selectWindow(name);
		run("Duplicate...", "title=contour duplicate channels="+contour_chan);
		
		
		cell_added=false;
		setTool(3);
		
		if (File.exists(dir_output+subname+"_CellROIs.zip")) {
			roiManager("open", dir_output+subname+"_CellROIs.zip");
			nCells=roiManager("count");
			selectWindow("contour");
			waitForUser("Draw cells and add to ROI Manager");
			if (roiManager("count")>nCells) {
				cell_added=true;
			}
		}
		
		else {
			while (roiManager("count")<1) {
				selectWindow("contour");
				waitForUser("Draw cell and add to ROI Manager");
			}
		}
			
		roiManager("deselect");
		roiManager("save", dir_output+subname+"_CellROIs.zip");
		
		
		selectWindow("Translocated_Prot");
		roiManager("deselect");
		roiManager("measure");
		
		NumbCell = (nResults); //The number of results = the number of cells detected
		PictureName=newArray(NumbCell);
		Area=newArray(NumbCell);//We create and fill the arrays with the information (name of image and cytosol intensity for each cell)
		TotMean=newArray(NumbCell);
		NucMean=newArray(NumbCell);

		for (cell_res = 0; cell_res < nResults; cell_res++) {
			PictureName[cell_res]=subname;
			Area[cell_res]=getResult("Area", cell_res);
			TotMean[cell_res]=getResult("Mean", cell_res);
		}	
				
		run("Clear Results");

		// Nuclei detection
		
		// If nucleus ROIs exists from previous analysis, and no cell has been added in this current analysis, open them
		if (File.exists(dir_output+subname+"_NucleiROIs.zip") && !cell_added) { 
			roiManager("reset");
			roiManager("open", dir_output+subname+"_NucleiROIs.zip");			
		}

		else {
			selectWindow("nucleus");
			run("Gaussian Blur...", "sigma="+sigma);
			setAutoThreshold("MaxEntropy dark");
			waitForUser("Adjust threshold for nuclei detection");
			
			for (nuc = 0; nuc < NumbCell; nuc++) {
				nuclei=roiManager("count");
				selectWindow("nucleus");
				roiManager("deselect");
				roiManager("select", nuc);
				run("Analyze Particles...", "size="+minSizeNuclei+"-Infinity add");
				while (roiManager("count")!=nuclei+1) {
					waitForUser("Check nucleus ROI");
				}

			}
	
			for (cell_del = 0; cell_del < NumbCell; cell_del++) {
				roiManager("select", 0);
				roiManager("delete");
			}
		}
		roiManager("deselect");
		roiManager("save", dir_output+subname+"_NucleiROIs.zip");
		
		selectWindow("Translocated_Prot");
		roiManager("deselect");
		roiManager("measure");
			
		
		NumbNuclei=nResults;
		
		for (nuc_res=0; nuc_res<NumbNuclei; nuc_res++) {
			NucMean[nuc_res] = getResult("Mean", nuc_res);
		}

		close("*");		
		run("Clear Results");
		roiManager("reset");		
	}
	
	else {
		continue
	}

	//We add the informations obtained for this image to the master arrays that will recapitulate all the informations for the dir_input 
	Name=Array.concat(Name,PictureName);
	CellArea=Array.concat(CellArea,Area);
	NucleiMean=Array.concat(NucleiMean,NucMean);
	TotalMean=Array.concat(TotalMean,TotMean);	
	
	run("Clear Results");
	roiManager("reset");		
}
	

if (nResults!=0) {
	run("Clear Results");	
}	


//We create a result table with all the informations gathered during the dir_input analysis and we'll save it as a master result Excel file
for (j = 0; j < lengthOf(Name); j++) {
	setResult("Name", j, Name[j]);
	setResult("Cell area", j, CellArea[j]);
	setResult("Nuclei Mean Intensity", j, NucleiMean[j]);
	setResult("Cytoplasm Mean Intensity", j, TotalMean[j]);
	setResult("Ratio Mean Intensity Nuc/Total", j, NucleiMean[j]/TotalMean[j]);	
}	

saveAs("Results", dir_output+condition+"_TranslocatedProtein_IntensityRatio.xls");
selectWindow("Log");
saveAs("text", dir_output+condition+"_Z-Projections.txt");
run("Close All");


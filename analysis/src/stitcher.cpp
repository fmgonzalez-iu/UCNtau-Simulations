#include "../inc/stitcher.hpp"

std::vector<noabsResult> stitch_noabs_to_noabsClean(std::vector<noabsCleanDagResult> evts_in) {
	// Convert noabsCleanDagResult to noabsResult
	
	std::vector<noabsResult> evts_out;
	evts_out.reserve(evts_in.size()); // Need to properly allocate size here
	for (auto event = evts_in.begin(); event<evts_in.end(); event++) {
		
		noabsResult tmp_evt;
		tmp_evt.energy = (*event).energy;
		tmp_evt.theta  = (*event).theta;
				
		for (size_t i = 0; i < NRECORDS; i++) {
			tmp_evt.times[i]  = (*event).times[i];
			tmp_evt.ePerps[i] = (*event).ePerps[i];
			tmp_evt.zetas[i]  = 0.0; // This is the patch here
		}
		evts_out.push_back(tmp_evt);
		
	}
	printf("Sucessfully stitched file!\n");
	return evts_out;
}

TH1D* stitch_weighted_vec_to_root(std::vector<weightedBin> evts_in) {
	
	char hName[24]; // This is our ROOT histogram
    sprintf(hName, "MCEvents");
    TH1D* hist = new TH1D(hName, "Monte Carlo", (int)evts_in.size(), 0, (int)evts_in.size());
    
    double dipEnds[NDIPS] = DIPTIMES; // Dip times
    printf("dipEnds loaded\n");
	
	for (size_t i = 0; i < evts_in.size(); i++) {
		hist->SetBinContent(i, evts_in[i].wgt);// Fill ROOT hist
	}
	return hist;
	
}

TH1D* stitch_fixed_result_to_root(std::vector<fixedResult> evts_in) {
	
	char hName[24]; // this is our ROOT histogram
	sprintf(hName, "MCEvents");
	TH1D* hist = new TH1D(hName, "Monte Carlo", 310 , 0, 310);
	
	for (size_t i = 0; i < evts_in.size(); i++) {
		//printf("%f\n",evts_in[i].t);
		if(evts_in[i].zOff == -5 && evts_in[i].nHit < 1) { // Put conditions here
		//if(evts_in[i].zOff == -5) {
			
			hist->Fill(evts_in[i].t);
		}
	}
	return hist;
}

TH2D* create_root_2D(std::vector<std::vector<double>> state) {
	// 2D slice histogram
	char hName[24]; // Name the histogram
	sprintf(hName, "MCPositions");
	TH2D* hist = new TH2D(hName, "Monte Carlo", 200,-1.0,1.0,200,-1.0,1.0); // 1cm bins, 
	
	for (size_t i = 0; i < state.size(); i++) {
		if (state[i][0] > 0) { // Put some "time" condition here.
			hist->Fill(state[i][1],state[i][2]); // fill histogram
		}
	}
	return hist;
}
	

TH3D* create_root_3D(std::vector<std::vector<double>> state) {
	
	char hName[24]; // Name the histogram
	sprintf(hName, "MCPositions");
	TH3D* hist = new TH3D(hName, "Monte Carlo", 200,-1.0,1.0,200,-1.0,1.0,50,-1.5,-1.0);
	
	for (size_t i = 0; i < state.size(); i++) {
		if (state[i][0] > 0) { // Put some "time" condition here.
		hist->Fill(state[i][1],state[i][2],state[i][3]); // Fill histogram
		}
	}
	return hist;
	
}

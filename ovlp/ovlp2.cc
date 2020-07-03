#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>

#include <unordered_map>
#include <string.h>

using namespace std;

#include "./../om_set1/msfl.cpp"
#include "./../om_set1/m_read.cpp"
#include "./../om_set1/scoring.cpp"
#include "./../om_set1/alignment.cpp"


int main(int argc, char *argv[])
{
  if(argc != 5){
    cerr<<"Usage:"<<endl;
    cerr<<"ovlp2 <Map_File> <Candidate_File> <Output> <Detailed_Output>"<<endl;
  }
  else{
    ofstream det_str(argv[4]);
    assert(det_str.good());

    //ifstream ifs("./../datasets/pestis/rough_map_2");
    //ifstream ifs("./../datasets/pestis/rough_maps.v2");
    //ifstream ifs("./../datasets/pestis/YPestis.kim.rfmap");
    //ifstream ifs("./../datasets/pestis/rm2.cons");
    //ifstream ifs("./../datasets/Ypestis_XhoI.maps.ex"); //file with the om data
    //ifstream ifs("./../datasets/pestis/Ypestis_XhoI.maps_original");
    //ifstream ifs("./../datasets/susan/chr19/mole_candidate_cut5_chr19.maps");
    //ifstream ifs("./../datasets/pestis/Ypestisoptical.maps");
    //ifstream ifs("./../datasets/Ypestis_XhoIcontigged.maps.short");
    //ifstream ifs("./../datasets/Ypestis_XhoIcontigged.maps");
    //ifstream ifs("./datasets/Ypestis.perm");
    //ifstream ifs("./datasets/Ypestis_XhoIcontigged.maps.77.nonflip");
    //ifstream ifs("./datasets/ch13.om.sim");
    //ifstream ifs("./../datasets/pfluor/Pfluorescensoptical.cut");
    ifstream ifs;
    ifs.open(argv[1]);
    if(!ifs.good()){
      cerr<<"wrong maps file name: "<<argv[1]<<endl;
      assert(false);
    }
    //assert(ifs.good());
    
    om_read_collection maps(ifs);

    ifs.close();

    unordered_map<string, int> name2id;
    for(int i=0; i<maps.collection.size(); i++){
      name2id[maps.collection[i].read_name] = i;
    }
    
    //ifstream ref_if("./../datasets/pestis/YPestis.kim.rfmap.ext");
    //ifstream ref_if("./../datasets/susan/chr19/chr19_silico.maps");
    //ifstream ref_if("./datasets/YPestis.kim.rfmap");
    //ifstream ref_if("./../datasets/YPestis.kim.rfmap.ex");
    //ifstream ref_if("./datasets/ch13.ref");
    //assert(ref_if.good());

    remove(argv[3]);

    ofstream times;
    times.open("valuev_row_times.csv");
    ofstream ovlps;
    ovlps.open(argv[3]);    

    ifstream pair_ifs;
    pair_ifs.open(argv[2]);
    if(!pair_ifs.good()){
      cerr<<"wrong pair file name: "<<argv[1]<<endl;
      assert(false);
    }

    // unordered_map<int,string> pairid2name;
    // string line;
    // int j = 0;
    // while(getline(pair_ifs, line)) {
    //   char *sl = new char[256];
    //   strncpy(sl, line.c_str(), 256);
    //   sl[255] = '\0';
    //   for(int i = 0; i < strlen(sl); i++) {
    // 	if (sl[i] == ':') {
    // 	  sl[i] = '\0';
    // 	  break;
    // 	}
    //   }
    //   pairid2name[j] = (string) sl;
    //   j++;
    // }

    // pair_ifs.clear();
    // pair_ifs.seekg(0);

    scoring_params sp(.2,1.2,.9,3,17.43,0.58, 0.0015, 0.8, 1, 3);
    //sp.init();

    int stored = 0;
    string line;
    while(getline(pair_ifs, line)) {
      clock_t start_time = clock();
      char name[256];
      const char *sl = line.c_str();

      int i = 0;
      for(; i < strlen(sl); i++) {
	if (sl[i] == ':') {
	  name[i] = '\0';
	  break;
	}
	name[i] = sl[i];
      }

      i+=2;
      int id1;
      int id2;
      char ori;

      if (name2id.count(name) != 1) {
	cerr << "Invalid Rmap name " << name << endl;
	assert(false);
      }
      
      id1 = name2id[name];
      while(i < strlen(sl)) {
	char name2[256];
	int j = 0;
	for(;i < strlen(sl); i++, j++) {
	  if (sl[i] == ',') {
	    name2[j] = '\0';
	    break;
	  }
	  name2[j] = sl[i];
	}
	if (name2id.count(name2) != 1) {
	  cerr << "Invalid Rmap name " << name2 << endl;
	  assert(false);
	}
	id2 = name2id[name2];
	while(i < strlen(sl) && sl[i] != ',') i++;
	i++;
	ori = sl[i];
	i++;
	if (sl[i] != ' ' && sl[i] != '\0' && sl[i] != '\n') {
	  cerr<<"malformed line: "<<line<< " " <<  i << endl;
	  assert(false);
	}
	i++;
	
	if(id1!=id2 && maps.collection[id1].map_read.size() > 5
	   && maps.collection[id2].map_read.size() > 5){

	  om_read tar_map = maps.collection[id1];
	  om_read for_map = maps.collection[id2];
	  if (ori == 'r') {
	    for_map = for_map.reverse();
	  }
	  
	  rm_alignment for_alignment(tar_map, for_map, sp);

	  //for_alignment.localized_overlap_alignment
	  //  (0,tar_map.map_read.size(),0,for_map.map_read.size());

	  for_alignment.optimized_overlap_alignment();

	  //for_alignment.overlap_alignment();

	  for_alignment.overlap_t_score();
	      
	  double for_score = for_alignment.Smax;

	  double for_t_score = for_alignment.Tmax;

	  //double for_p_value = for_alignment.ovlp_p_value();

	  double for_ovlp_size = for_alignment.ovlp_size();

	  // cout<<"fs: "<<for_score<<" ft: "<<for_t_score;
	  // //cout<<" p_v: "<<for_p_value;
	  // cout<<endl;

	  double score_thresh = 21; // originally 25
	  double t_score_thresh = 7; // originally 8
	  double t_mult = 0;


	  if(for_t_score > t_score_thresh &&
	     for_score > score_thresh){
	    //&& for_t_score > t_mult*for_ovlp_size){
	    stored++;
	    int ovlp_start1 = 
	      for_alignment.ref_restr_al_sites
	      [for_alignment.ref_restr_al_sites.size()-1];
	    int ovlp_end1 = for_alignment.ref_restr_al_sites[0];

	    int ovlp_start2 = 
	      for_alignment.tar_restr_al_sites
	      [for_alignment.tar_restr_al_sites.size()-1];
	    int ovlp_end2 = for_alignment.tar_restr_al_sites[0];

	    ovlps<<tar_map.read_name.c_str();	    
	    ovlps<<" "<<for_map.read_name.c_str();
	    ovlps<<" "<<tar_map.map_read.size();
	    ovlps<<" "<<for_map.map_read.size();
	    ovlps<<" 1 1 "<<for_score;
	    ovlps<<" "<<for_t_score<<endl;//" ";
	    for(int k=for_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
	      if(k!=for_alignment.ref_restr_al_sites.size()-1)
		ovlps<<" ";
	      ovlps<<for_alignment.ref_restr_al_sites[k];
	      ovlps<<" ";
	      ovlps<<for_alignment.tar_restr_al_sites[k];
	    }
	    //ovlps<<tar_map.map_read.size()<<" ";
	    //ovlps<<for_map.map_read.size()<<" ";
	    //ovlps<<ovlp_start1<<" "<<ovlp_end1<<" ";
	    //ovlps<<ovlp_start2<<" "<<ovlp_end2<<endl;
	    ovlps<<endl<<endl;

	    //if(for_score>score_thresh)
	    //if(for_t_score > t_score_thresh)
	    for_alignment.output_alignment(cout);
	    for_alignment.output_alignment(det_str);
	  }
	}
	  
	if (sl[i] == '\n')
	  i++;
	  
      }

      clock_t end_time = clock();
      times << i << ", " << end_time - start_time  << std::endl;
    }
    
    ovlps.close();
    det_str.close();
    times.close();
    
  }
  return 0;  
}
  

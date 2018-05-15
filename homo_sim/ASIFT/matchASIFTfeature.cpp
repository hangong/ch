// Copyright (c) 2008-2011, Guoshen Yu <yu@cmap.polytechnique.fr>
// Copyright (c) 2008-2011, Jean-Michel Morel <morel@cmla.ens-cachan.fr>
//
// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// Jean-Michel Morel and Guoshen Yu, Method and device for the invariant 
// affine recognition recognition of shapes (WO/2009/150361), patent pending. 
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.
//
// 
//*----------------------------- demo_ASIFT  --------------------------------*/
// Detect corresponding points in two images with the ASIFT method. 

// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
// 
// Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
//            Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 
// Reference: ASIFT online demo (You can try ASIFT with your own images online.) 
//			  http://www.ipol.im/pub/algo/my_affine_sift/
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "demo_lib_sift.h"
#include "io_png/io_png.h"

#include "library.h"
#include "frot.h"
#include "fproj.h"
#include "compute_asift_keypoints.h"
#include "compute_asift_matches.h"

# define IM_X 800
# define IM_Y 600

void initkey(int num_tilt, vector< vector<keypointslist> >& keys_all);

int main(int argc, char **argv)
{			
	
    if ((argc != 4)) {
        std::cerr << " ******************************************************************************* " << std::endl
				  << " ***************************  ASIFT feature matching  **************************** " << std::endl
				  << " ******************************************************************************* " << std::endl
				  << "Usage: " << argv[0] << " key1.txt key2.txt " << std::endl
										  << " matchings.txt " << std::endl
										  << "- matchings.txt: coordinates of matched points (col1, row1, col2, row2). " << std::endl
										  << "- keys1.txt keys2.txt: ASIFT keypoints of the two images." << std::endl
   				  << " ******************************************************************************* " << std::endl
				  << " *********************  Han Gong, 2015 ******************** " << std::endl
				  << " *********************  Jean-Michel Morel, Guoshen Yu, 2010 ******************** " << std::endl
				  << " ******************************************************************************* " << std::endl;
        return 1;
    }
	
	//////////////////////////////////////////////// Input
	int wS1, hS1, wS2, hS2;

	///// Compute ASIFT keypoints
	// number N of tilts to simulate t = 1, \sqrt{2}, (\sqrt{2})^2, ..., {\sqrt{2}}^(N-1)
	int num_of_tilts1 = 7;
	int num_of_tilts2 = 7;
	int verb = 0;
    int rnum;
    int num_keys1, num_keys2;
	// Define the SIFT parameters
	siftPar siftparameters;	
	default_sift_parameters(siftparameters);

	vector< vector< keypointslist > > keys1;		
	vector< vector< keypointslist > > keys2;	
	
	// Write all the keypoints (row, col, scale, orientation, desciptor (128 integers)) to 
	// the file argv[6] (so that the users can match the keypoints with their own matching algorithm if they wish to)
	// keypoints in the 1st image
    initkey(num_of_tilts1, keys1);
	std::ifstream file_key1(argv[1]);
	if (file_key1.is_open())
	{
		// Follow the same convention of David Lowe: 
		// the first line contains the number of keypoints and the length of the desciptors (128)
		file_key1 >> num_keys1 >> wS1 >> hS1;
		for (int tt = 0; tt < (int) keys1.size(); tt++)
		{
			for (int rr = 0; rr < (int) keys1[tt].size(); rr++)
			{
                file_key1 >> rnum;
                keys1[tt][rr] = std::vector < keypoint > (rnum);
				keypointslist::iterator ptr = keys1[tt][rr].begin();
				for(int i=0; i < (int) keys1[tt][rr].size(); i++, ptr++)	
				{
					file_key1 >> ptr->x >> ptr->y >> ptr->scale >> ptr->angle;
					
					for (int ii = 0; ii < (int) VecLength; ii++)
					{
						file_key1 >> ptr->vec[ii];
					}
				}
			}	
		}
	}
	else 
	{
		std::cerr << "Unable to open the file keys1."; 
	}

	file_key1.close();
	
	// keypoints in the 2nd image
    initkey(num_of_tilts2, keys2);
	std::ifstream file_key2(argv[2]);
	if (file_key2.is_open())
	{
		// Follow the same convention of David Lowe: 
		// the first line contains the number of keypoints and the length of the desciptors (128)
		file_key2 >> num_keys2 >> wS2 >> hS2;
		for (int tt = 0; tt < (int) keys2.size(); tt++)
		{
			for (int rr = 0; rr < (int) keys2[tt].size(); rr++)
			{
                file_key2 >> rnum;
                keys2[tt][rr] = std::vector < keypoint > (rnum);
				keypointslist::iterator ptr = keys2[tt][rr].begin();
				for(int i=0; i < (int) keys2[tt][rr].size(); i++, ptr++)	
				{
					file_key2 >> ptr->x >> ptr->y >> ptr->scale >> ptr->angle;
					
					for (int ii = 0; ii < (int) VecLength; ii++)
					{
						file_key2 >> ptr->vec[ii];
					}
				}
			}	
		}
	}
	else 
	{
		std::cerr << "Unable to open the file keys1."; 
	}

	file_key2.close();

	//// Match ASIFT keypoints
	time_t tstart, tend;	
	int num_matchings;
	matchingslist matchings;	
	cout << "Matching the keypoints..." << endl;
	tstart = time(0);
	num_matchings = compute_asift_matches(num_of_tilts1, num_of_tilts2, wS1, hS1, wS2, 
										  hS2, verb, keys1, keys2, matchings, siftparameters);
	tend = time(0);
	cout << "Keypoints matching accomplished in " << difftime(tend, tstart) << " seconds." << endl;

	////// Write the coordinates of the matched points (row1, col1, row2, col2) to the file argv[5]
	std::ofstream file(argv[3]);
	if (file.is_open())
	{		
		// Write the number of matchings in the first line
		file << num_matchings << std::endl;
		
		matchingslist::iterator ptr = matchings.begin();
		for(int i=0; i < (int) matchings.size(); i++, ptr++)		
		{
			file << ptr->first.x << "  " << ptr->first.y << "  " <<  ptr->second.x << 
			"  " <<  ptr->second.y << std::endl;
		}		
	}
	else 
	{
		std::cerr << "Unable to open the file matchings."; 
	}

	file.close();

    return 0;
}

void initkey(int num_tilt, vector< vector<keypointslist> >& keys_all)
{

  float t_min, t_k;
  int tt, num_rot_t2;
  int counter_sim=0;

  num_rot_t2 = 10;

  t_min = 1;
  t_k = sqrt(2.);

  keys_all = std::vector< vector< keypointslist > >(num_tilt);	

  for (tt = 1; tt <= num_tilt; tt++)
  {
    float t = t_min * pow(t_k, tt-1);

    if ( t == 1 )
    {
      counter_sim ++;

      keys_all[tt-1] = std::vector< keypointslist >(1);
    }
    else
    {
      int num_rot1 = round(num_rot_t2*t/2);        
      if ( num_rot1%2 == 1 )
      {
        num_rot1 = num_rot1 + 1;
      }
      num_rot1 = num_rot1 / 2;
      counter_sim +=  num_rot1;

      keys_all[tt-1] = std::vector< keypointslist >(num_rot1);	
    }         		
  }

}


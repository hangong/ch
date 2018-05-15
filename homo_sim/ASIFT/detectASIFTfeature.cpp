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

# define IM_X 800
# define IM_Y 600

int main(int argc, char **argv)
{			
	
    if ((argc != 3) && (argc != 4)) {
        std::cerr << " ******************************************************************************* " << std::endl
				  << " ***************************  ASIFT image feature  **************************** " << std::endl
				  << " ******************************************************************************* " << std::endl
				  << "Usage: " << argv[0] << " imgIn1.png  " << std::endl
										  << " keys1.txt [Resize option: 0/1] " << std::endl
									      << "- imgIn1.png: input images (in PNG format). " << std::endl
										  << "- keys1.txt: ASIFT keypoints of the two images." << std::endl
										  << "- [optional 0/1]. 1: input images resize to 800x600 (default). 0: no resize. " << std::endl 
   				  << " ******************************************************************************* " << std::endl
				  << " *********************  Han Gong, 2015 ******************** " << std::endl
				  << " *********************  Jean-Michel Morel, Guoshen Yu, 2010 ******************** " << std::endl
				  << " ******************************************************************************* " << std::endl;
        return 1;
    }
	
	//////////////////////////////////////////////// Input
	// Read image1
    float * iarr1;
    size_t w1, h1;
    if (NULL == (iarr1 = read_png_f32_gray(argv[1], &w1, &h1))) {
        std::cerr << "Unable to load image file " << argv[1] << std::endl;
        return 1;
    }
    std::vector<float> ipixels1(iarr1, iarr1 + w1 * h1);
	free(iarr1); /*memcheck*/
	
	///// Resize the images to area wS*hW in remaining the apsect-ratio	
	///// Resize if the resize flag is not set or if the flag is set unequal to 0
	float wS = IM_X;
	float hS = IM_Y;
	
	float zoom1=0;	
	int wS1=0, hS1=0;
	vector<float> ipixels1_zoom;	
		
	int flag_resize = 1;
	if (argc == 4)
	{	
		flag_resize = atoi(argv[3]);
	}
	
	if ((argc == 3) || (flag_resize != 0))
	{
		cout << "WARNING: The input images are resized to " << wS << "x" << hS << " for ASIFT. " << endl 
		<< "         But the results will be normalized to the original image size." << endl << endl;
		
		float InitSigma_aa = 1.6;
		
		float fproj_p, fproj_bg;
		char fproj_i;
		float *fproj_x4, *fproj_y4;
		int fproj_o;

		fproj_o = 3;
		fproj_p = 0;
		fproj_i = 0;
		fproj_bg = 0;
		fproj_x4 = 0;
		fproj_y4 = 0;
				
		float areaS = wS * hS;

		// Resize image 1 
		float area1 = w1 * h1;
		zoom1 = sqrt(area1/areaS);
		
		wS1 = (int) (w1 / zoom1);
		hS1 = (int) (h1 / zoom1);
		
		int fproj_sx = wS1;
		int fproj_sy = hS1;     
		
		float fproj_x1 = 0;
		float fproj_y1 = 0;
		float fproj_x2 = wS1;
		float fproj_y2 = 0;
		float fproj_x3 = 0;	     
		float fproj_y3 = hS1;
		
		/* Anti-aliasing filtering along vertical direction */
		if ( zoom1 > 1 )
		{
			float sigma_aa = InitSigma_aa * zoom1 / 2;
			GaussianBlur1D(ipixels1,w1,h1,sigma_aa,1);
			GaussianBlur1D(ipixels1,w1,h1,sigma_aa,0);
		}
			
		// simulate a tilt: subsample the image along the vertical axis by a factor of t.
		ipixels1_zoom.resize(wS1*hS1);
		fproj (ipixels1, ipixels1_zoom, w1, h1, &fproj_sx, &fproj_sy, &fproj_bg, &fproj_o, &fproj_p, 
			   &fproj_i , fproj_x1 , fproj_y1 , fproj_x2 , fproj_y2 , fproj_x3 , fproj_y3, fproj_x4, fproj_y4); 
	}
	else 
	{
		ipixels1_zoom.resize(w1*h1);	
		ipixels1_zoom = ipixels1;
		wS1 = w1;
		hS1 = h1;
		zoom1 = 1;
	}

	
	///// Compute ASIFT keypoints
	// number N of tilts to simulate t = 1, \sqrt{2}, (\sqrt{2})^2, ..., {\sqrt{2}}^(N-1)
	int num_of_tilts1 = 7;
	int verb = 0;
	// Define the SIFT parameters
	siftPar siftparameters;	
	default_sift_parameters(siftparameters);

	vector< vector< keypointslist > > keys1;		
	
	int num_keys1;
	
	cout << "Computing keypoints on the image..." << endl;
	time_t tstart, tend;	
	tstart = time(0);

	num_keys1 = compute_asift_keypoints(ipixels1_zoom, wS1, hS1, num_of_tilts1, verb, keys1, siftparameters);
	
	tend = time(0);
	cout << "Keypoints computation accomplished in " << difftime(tend, tstart) << " seconds." << endl;
	
	// Write all the keypoints (row, col, scale, orientation, desciptor (128 integers)) to 
	// the file argv[6] (so that the users can match the keypoints with their own matching algorithm if they wish to)
	// keypoints in the 1st image
	std::ofstream file_key1(argv[2]);
	if (file_key1.is_open())
	{
		file_key1 << num_keys1 << " " << w1 << "  " << h1 << "  " << std::endl;
		// Follow the same convention of David Lowe: 
		// the first line contains the number of keypoints and the length of the desciptors (128)
		for (int tt = 0; tt < (int) keys1.size(); tt++)
		{
			for (int rr = 0; rr < (int) keys1[tt].size(); rr++)
			{
				keypointslist::iterator ptr = keys1[tt][rr].begin();
                file_key1 << keys1[tt][rr].size() <<endl;
				for(int i=0; i < (int) keys1[tt][rr].size(); i++, ptr++)	
				{
					file_key1 << zoom1*ptr->x << "  " << zoom1*ptr->y << "  " << zoom1*ptr->scale << "  " << ptr->angle;
					
					for (int ii = 0; ii < (int) VecLength; ii++)
					{
						file_key1 << "  " << ptr->vec[ii];
					}
					
					file_key1 << std::endl;
				}
			}	
		}
	}
	else 
	{
		std::cerr << "Unable to open the file keys."; 
	}

	file_key1.close();
	
    return 0;
}

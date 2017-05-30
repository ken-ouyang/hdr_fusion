#include <iostream>
#include <string>
#include "hdrFusion.hpp"
#include <sstream>
#include <cv.h>
#include <highgui.h>


using namespace std;

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

//input path,output path, saturation, exposure, contrast
int main(int argc, char* argv[])
{
    Mat outImage;
    if(argc < 6){
      return 0;
    }
    int n=argc-5;
//    std::cout << argc << std::endl;
    // for (int i = 0; i < argc; i++) {
        // std::cout << argv[i] << std::endl;
    // }
    vector<Mat> images(n);
    for(int i=0; i<n; i++){
       images[i] = imread(argv[i+1],1);
       //images[i] = imread("/home/ken/Documents/Test/HDRFusion/TestImg_Align/1_1.jpg", 1);
    }
    double saturation = atof(argv[n+2]);
    double exposure = atof(argv[n+3]);
    double contrast = atof(argv[n+4]);
    // saturation, exposure , contrast
    clock_t start_timer, end_timer;

    start_timer = clock();
    hdrFusion(images, outImage, saturation, exposure, contrast);
    end_timer = clock();

    std::cout << "hdrFusion Cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC * 1000 << " ms"<< std::endl;

    imwrite(argv[n+1], outImage);

    return 0;
}

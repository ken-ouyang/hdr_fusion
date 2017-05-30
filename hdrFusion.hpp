
//-------------------------------------------
//
// HDR_Fusion
//
// Created by OUYANG Hao @2016/6/8
//
//-------------------------------------------
#ifndef _hdrFusion_HPP_
#define _hdrFusion_HPP_

#include <opencv2/core/core.hpp>

using namespace cv;



int qualityEstimation(Mat &input_image, Mat &weight_map, double weight_saturation, double weight_exposure, double weight_contrast);

//HDR Fusion
//input_image, a set of image with different exposure
//output_image, hdr Results
//weight_saturation, the weight of saturation in qualityEstimation process
//weight_exposure, the weight of exposure in qualityEstimation process
//weight_contrast, the weight of contrast in qualityEstimation process

int hdrFusion(vector<Mat> &input_image, Mat &output_image, double weight_saturation, double weight_exposure, double weight_contrast);


#endif

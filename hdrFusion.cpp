
//-------------------------------------------
// HDR Fusion
//
// Created by OUYANG Hao @2016/6/7
//
//-------------------------------------------

#include "hdrFusion.hpp"
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "guided_filter.hpp"
#define TIMING
using namespace cv;

//int x=0;
//Quality Estimate For Fusion
int qualityEstimation(Mat &input_image, Mat &weight_map, double weight_saturation, double weight_exposure, double weight_contrast)
{
    CV_Assert(input_image.data != NULL);
    CV_Assert(input_image.channels()==1 || input_image.channels()==3 || input_image.channels()==4);
    Mat input;
    input_image.convertTo(input, CV_32F, 1.0/255.0);
    int width = input_image.cols;
    int height = input_image.rows;
    vector<Mat> channels(3);
    split(input, channels);
    //x=x+1;

    //Caculate Saturation
    Mat saturation_map;
    Mat mean_map;
    Mat temp1, temp2, temp3;

    if (input_image.channels()==1)
    {
        saturation_map = Mat::ones(height, width, CV_32F);
        weight_saturation = 0;
        mean_map = input;
    }
    else
    {
        mean_map = (channels[0] + channels[1] + channels[2]) / 3;
        pow((channels[0] - mean_map), 2, temp1);
        pow((channels[1] - mean_map), 2, temp2);
        pow((channels[2] - mean_map), 2, temp3);
        pow((temp1 + temp2 + temp3) / 3, 0.5 ,saturation_map);
    }

    //Caculate Exposure

#ifdef TIMING
    clock_t start_timer, end_timer;
    start_timer = clock();
#endif

    float exposure = 0.5;
    Mat exposure_map;
    vector<Mat> temp_channels(3);
    double std = 0.2;
    if (input_image.channels() == 1)
    {
        pow((channels[0] - exposure) , 2, temp1);
        exp(-temp1 / (2 * std * std), exposure_map);
    }
    else
    {
        pow((channels[0] - exposure) , 2, temp1);
        exp(-temp1 / (2 * std * std), temp_channels[0]);
        pow((channels[1] - exposure) , 2, temp2);
        exp(-temp2 / (2 * std * std), temp_channels[1]);
        pow((channels[2] - exposure) , 2, temp3);
        exp(-temp3 / (2 * std * std), temp_channels[2]);
        exposure_map = (temp_channels[0].mul(temp_channels[1])).mul(temp_channels[2]);
        // pow((temp_channels[0].mul(temp_channels[1])).mul(temp_channels[2]), 1.0/3.0, exposure_map);
    }

#ifdef TIMING
    end_timer = clock();
    std::cout << "exposure estimation cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC * 1000 << " ms" <<  std::endl;
#endif
    //Caculate Contrast Using Laplacian Filter

    Mat contrast_map;
    Laplacian(mean_map * 255, contrast_map, CV_32F, 3, 1, 0, BORDER_DEFAULT);
    convertScaleAbs(contrast_map, contrast_map);
    contrast_map.convertTo(contrast_map, CV_32F, 1.0/255.0);

    //Caculate the overall grade
    pow(saturation_map, weight_saturation, saturation_map);
    pow(exposure_map, weight_exposure, exposure_map);
    pow(contrast_map, weight_contrast, contrast_map);

    imwrite("./TestImg_Align/con_map.jpg",contrast_map*255);
    // imwrite("./TestImg_Align/sat_map"+patch::to_string(x)+".jpg",saturation_map*255);
    imwrite("./TestImg_Align/exp_map.jpg",exposure_map*255);

    double x = weight_contrast + weight_saturation + weight_exposure;
    pow((saturation_map.mul(exposure_map)).mul(contrast_map), 1.0/x, weight_map);
    // weight_map = (exposure_map.mul(contrast_map)).mul(saturation_map);

    imwrite("./TestImg_Align/weight_map.jpg",weight_map*255);


}


void weight_adjustment_by_ZNCC(vector<Mat> & input_image, vector<Mat> & weight_map, int weight_ZNCC)
{
    // convert to YUV channels
    int image_num = input_image.size();
    int ref_num = image_num / 2;
    vector<Mat> y_channel_set;
    for (int i = 0; i < image_num; i ++)
    {
        Mat temp;
        vector<Mat> temp_channels;
        cvtColor(input_image[i], temp, CV_BGR2YCrCb);
        split(temp, temp_channels);
        y_channel_set.push_back(temp_channels[0]);
    }


    // //////////////////////////////////////////////////////////////////

    // reject the patches which do not satisfy the brightness requirement

    // //////////////////////////////////////////////////////////////////

    int patch_size_c = 11;
    Mat temp_diff;
    Mat comp_map;
    bool flag = true;
    // std::cout << "width "<< y_channel_set[0].cols << std::endl;
    // std::cout << "height "<< y_channel_set[0].rows << std::endl;


    for (int r = 0; r < input_image[0].rows - patch_size_c; r = r + patch_size_c)
    {
        for (int c = 0; c < input_image[0].cols - patch_size_c; c = c + patch_size_c)
        {
            // std::cout << "r " << r << " c " << c << std::endl;
            // Mat roi_r(y_channel_set[ref_num], Rect(r, c, patch_size_c, patch_size_c));
            for (int i = 0; i < ref_num; i ++)
            {
                // Mat roi_c(y_channel_set[i], Rect(r, c, patch_size_c, patch_size_c));
                // temp_diff = roi_c - roi_r;
                // compare(temp_diff, 0, comp_map, CMP_LT);
                flag = true;
                for (int r1 = r; r1 < r + patch_size_c; r1 ++)
                {
                    if (flag == false)
                    {
                        break;
                    }
                    for (int c1 = c; c1 < c + patch_size_c; c1 ++)
                    {
                        if (y_channel_set[i].at<float>(r1, c1) - y_channel_set[ref_num].at<float>(r1, c1) > 0)
                        {
                            // Mat roi_w(weight_map[i], Rect(r, c, patch_size_c, patch_size_c)) = Scalar(0);
                            // std::cout << "width " << weight_map[i].cols << std::endl;
                            // std::cout << "height " << weight_map[i].rows << std::endl;
                            // std::cout << "r " << r << " c " << c << std::endl;
                            Mat roi_w(weight_map[i], Rect(c, r, patch_size_c, patch_size_c));
                            // std::cout << "testing " << std::endl;
                            roi_w = Scalar(0);
                            flag = false;
                            break;
                        }
                    }
                }


            }

            for (int i = ref_num + 1; i < image_num; i ++)
            {
                // Mat roi_c(y_channel_set[i], Rect(r, c, patch_size_c, patch_size_c));
                // temp_diff = roi_c - roi_r;
                // compare(temp_diff, 0, comp_map, CMP_LT);
                flag = true;
                for (int r1 = r; r1 < r + patch_size_c; r1 ++)
                {
                    if (flag == false)
                    {
                        break;
                    }
                    for (int c1 = c; c1 < c + patch_size_c; c1 ++)
                    {
                        if (y_channel_set[i].at<float>(r1, c1) - y_channel_set[ref_num].at<float>(r1, c1) < 0)
                        {
                            Mat roi_w(weight_map[i], Rect(c, r, patch_size_c, patch_size_c));
                            roi_w = Scalar(0);
                            flag = false;
                            break;
                        }
                    }
                }

            }
        }
    }

    /////////////////////////////////////////////////////////////////////

    // adjust the weight by ZNCC

    /////////////////////////////////////////////////////////////////////

    int patch_size_z = 33;
    Mat diff_ref;
    Mat diff_i;
    float diff_ref_total;
    float diff_i_total;
    float diff_i_ref_total;
    Mat compare_tmp;
    for (int r = 0; r < input_image[0].rows - patch_size_z; r = r + patch_size_z)
    {
        for (int c = 0; c < input_image[0].cols - patch_size_z; c = c + patch_size_z)
        {
            // Calculate ZNCC for different patch

            Mat roi_r(y_channel_set[ref_num], Rect(c, r, patch_size_z, patch_size_z));
            float mean_ref = mean(roi_r)[0];
            diff_ref = roi_r - mean_ref;
            diff_ref_total = sum(diff_ref.mul(diff_ref))[0];

            // compare(abs(diff_ref), 0.01, compare_tmp, CMP_LE);
            // int sum_flat = sum(compare_tmp)[0] / 255;

            for (int i = 0; i < image_num; i ++)
            {
                if (i == ref_num)
                {
                    continue;
                }
                Mat roi_c(y_channel_set[i], Rect(c, r, patch_size_z, patch_size_z));
                float mean_cur = mean(roi_c)[0];
                diff_i = roi_c - mean_cur;
                diff_i_total = sum(diff_i.mul(diff_i))[0];
                diff_i_ref_total = sum(diff_i.mul(diff_ref))[0];
                float zncc;

                // handle extreme case when ref is under or over exposured

                // method 1 : Estimate the number of diff i

                // If the most part of the region is only of single color, then we can probobaly assume that this part do not have a ghost effect

                // int temp = countNonZero(diff_ref);


                if (diff_ref_total == 0 || diff_i_total == 0)
                // if (1)
                // if (sum_flat > patch_size_z * patch_size_z / 10 || diff_i_total == 0)
                {
                    zncc = 1;
                }
                else
                {
                    zncc = diff_i_ref_total / sqrt(diff_i_total * diff_ref_total) ;
                }
                // std::cout << "zncc " << zncc << std::endl;
                Mat roi_w(weight_map[i], Rect(c, r, patch_size_z, patch_size_z));
                roi_w = roi_w * pow(zncc, 30);
            }

        }
    }


}

// Note that all input_image are of the same size
// The smaller one of the height and weight is expected to be larger than 30 pixels
int hdrFusion(vector<Mat> &input_image, Mat &output_image, double weight_saturation, double weight_exposure, double weight_contrast)
{

    int image_num = input_image.size();
    int width = input_image[0].cols;
    int height = input_image[0].rows;
    vector<Mat> weight_map(image_num);
    int i, j;

#ifdef TIMING
    clock_t start_timer;
    clock_t end_timer;
#endif


#ifdef TIMING
    start_timer = clock();
#endif

    for (i=0; i<image_num; i++)
    {
        qualityEstimation(input_image[i], weight_map[i], weight_saturation, weight_exposure, weight_contrast);
        input_image[i].convertTo(input_image[i],CV_32F, 1.0/255.0);
    }

    // deghost method 1 : suppose the image set is ascending order of brightness
    int weight_ZNCC = 10;

    weight_adjustment_by_ZNCC(input_image, weight_map, weight_ZNCC);



#ifdef TIMING
    end_timer = clock();
    std::cout << " Quality Estimation Cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC * 1000 << " ms" <<  std::endl;
#endif


#ifdef TIMING
    start_timer = clock();
#endif

    int level = log2(min(width,height)/10);
    std::cout << "level " << level  << std::endl;
    if (level < 2)
    {
        std::cout << "Image Size is Too Small" << std::endl;
        return 0;
    }
    int pyr_height [level];
    int pyr_width [level];
    vector< vector<Mat> > weight_map_pyr;
    vector< vector<Mat> > input_image_pyr_gau;
    vector< vector<Mat> > input_image_pyr_lap;
    vector< vector<Mat> > input_image_pyr_temp;
    vector<Mat> output_image_pyr(level);
    weight_map_pyr.resize(image_num);
    input_image_pyr_temp.resize(image_num);
    input_image_pyr_lap.resize(image_num);
    input_image_pyr_gau.resize(image_num);
    for(i=0; i<image_num; i++)
    {
        weight_map_pyr[i].resize(level);
        input_image_pyr_gau[i].resize(level);
        input_image_pyr_lap[i].resize(level);
        input_image_pyr_temp[i].resize(level);
    }
    pyr_width[0] = width; pyr_height[0] = height;

#ifdef TIMING
    end_timer = clock();
    std::cout << " Pyramid Initialize Cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC * 1000 << " ms" <<  std::endl;
#endif

    //Construct Pyramid

#ifdef TIMING
    start_timer = clock();
#endif

    Mat weight_temp;
    Mat input_temp;
    for(j=0; j<image_num; j++)
    {
        //Construct Gaussian Pyramid for weight map
        weight_map_pyr[j][0] = weight_map[j];
        for(i=0; i<level-1; i++)
        {
            if(j==0)
            {
                pyr_width[i+1] = weight_map_pyr[j][i].cols / 2;
                pyr_height[i+1] = weight_map_pyr[j][i].rows / 2;
            }
            //Gaussian Pyramid
            //pyrDown(weight_map_pyr[j][i] , weight_map_pyr[j][i+1], Size(pyr_width[i+1], pyr_height[i+1]));
            //Bilateral Pyramid
            //bilateralFilter ( weight_map_pyr[j][i], weight_temp, 3, 5, 5);
            //guidedfilter
            fastguidedfilter(weight_map_pyr[j][i], weight_map_pyr[j][i], weight_temp, 3, 0.1, 5);
            resize(weight_temp, weight_map_pyr[j][i+1],  Size(pyr_width[i+1], pyr_height[i+1]), 0, 0, INTER_NEAREST);
            //imwrite("./TestImg_Align/pyr"+patch::to_string(i)+".jpg",weight_map_pyr[j][i] * 255);
        }

        //Construct Laplacian Pyramid for original image
        input_image_pyr_gau[j][0] = input_image[j] ;
        for(i=0; i<level-1; i++)
        {
            //Gaussian Pyramid
            //pyrDown(input_image_pyr_gau[j][i], input_image_pyr_gau[j][i+1], Size(pyr_width[i+1], pyr_height[i+1]));
            //Bilateral Pyramid
            //bilateralFilter(input_image_pyr_gau[j][i], input_image_pyr_temp[j][i], 3, 5, 5);
            //guidedfilter
            fastguidedfilterColor(input_image_pyr_gau[j][i], input_image_pyr_gau[j][i], input_image_pyr_temp[j][i], 3, 0.1, 5);
            resize(input_image_pyr_temp[j][i], input_image_pyr_gau[j][i+1], Size(pyr_width[i+1], pyr_height[i+1]), 0, 0, INTER_NEAREST);
            //pyrUp(input_image_pyr_gau[j][i+1], input_image_pyr_temp[j][i], Size(pyr_width[i], pyr_height[i]));
            resize(input_image_pyr_gau[j][i+1], input_image_pyr_temp[j][i], Size(pyr_width[i], pyr_height[i]), 0, 0, INTER_LINEAR);
            //bilateralFilter(input_image_pyr_temp[j][i], input_temp, 3, 5, 5);
            fastguidedfilterColor(input_image_pyr_temp[j][i], input_image_pyr_temp[j][i], input_temp, 3, 0.1, 5);
            if(i==level-2)
            {
                input_image_pyr_lap[j][i] = input_image_pyr_gau[j][i];
            }
            else
            {
                input_image_pyr_lap[j][i] = (input_image_pyr_gau[j][i] - input_temp);
                //input_image_pyr_lap[j][i] = (input_image_pyr_gau[j][i] - input_image_pyr_temp[j][i]);
            }
            //imwrite("./TestImg_Align/lapPyr"+patch::to_string(i)+patch::to_string(j) + ".jpg",(input_image_pyr_lap[j][i] + 0.5) * 255);
        }
        //imwrite("./TestImg_Align/lapPyr"+patch::to_string(level-2)+".jpg",(input_image_pyr_lap[j][level-2] ) * 255);
  }

#ifdef TIMING
  end_timer = clock();
  std::cout << " Construct Pyramid Cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC * 1000 << " ms"  <<  std::endl;
#endif

  //Combine the Pyramid

#ifdef TIMING
  start_timer = clock();
#endif

    Mat mask;
    vector<Mat> weight_sum(level);


    for (i=0; i<level-1; i++)
    {
        for (j=0; j<image_num; j++)
        {
            if(weight_sum[i].data != NULL)
            {
                weight_sum[i] = weight_sum[i] + weight_map_pyr[j][i];
            }
            else
            {
                weight_sum[i] = weight_map_pyr[j][i].clone();
            }
        }

        // handle the case for all weight is zero

        compare(weight_sum[i], 0, mask, CMP_EQ);
        weight_sum[i].setTo(1, mask);
        weight_map_pyr[0][i].setTo(1,mask);


        //imwrite("./TestImg_Align/weightPyr"+patch::to_string(i)+".jpg",weight_sum[i] );

        for (j=0; j<image_num; j++)
        {
            divide(weight_map_pyr[j][i],weight_sum[i],weight_map_pyr[j][i],1);
            //imwrite("./TestImg_Align/weightPyr2"+patch::to_string(i)+".jpg",weight_map_pyr[j][i] );
            cvtColor(weight_map_pyr[j][i],weight_map_pyr[j][i], CV_GRAY2BGR);
            //imwrite("./TestImg_Align/weightPyr2"+patch::to_string(j) + patch::to_string(i)+".jpg",weight_map_pyr[j][i]*255 );
            if (output_image_pyr[i].data != NULL)
            {
                output_image_pyr[i] = output_image_pyr[i] + weight_map_pyr[j][i].mul(input_image_pyr_lap[j][i]);
            }
            else
            {
                output_image_pyr[i] =  (weight_map_pyr[j][i]).mul(input_image_pyr_lap[j][i]);
            }
        }
        //imwrite("./TestImg_Align/outPyr"+patch::to_string(i)+".jpg",output_image_pyr[i] * 255);
    }

#ifdef TIMING
  end_timer = clock();
  std::cout << " Combine Pyramid Cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC  * 1000 << " ms" <<  std::endl;
#endif


    //Collapse the Pyramid

#ifdef TIMING
  start_timer = clock();
#endif
    Mat temp_mat = output_image_pyr[level - 2];
    Mat temp_mat2;
    for (i=0; i< level-2; i++)
    {
        //pyrUp(temp_mat, temp_mat, Size(pyr_width[level-3-i], pyr_height[level-3-i]));
        resize(temp_mat, temp_mat, Size(pyr_width[level-3-i], pyr_height[level-3-i]), 0, 0, INTER_LINEAR);
        fastguidedfilterColor(temp_mat, temp_mat, temp_mat2, 3, 0.1, 5);
        //bilateralFilter(temp_mat, temp_mat2, 3, 5, 5);
        temp_mat = temp_mat2 + output_image_pyr[level-3-i] ;
        //imwrite("./TestImg_Align/outPyr"+patch::to_string(i) + ".jpg", temp_mat*255);
    }
    output_image = temp_mat * 255;

    for (i=0; i<image_num; i++)
    {
        input_image[i].convertTo(input_image[i],CV_32F, 255.0);
    }

#ifdef TIMING
  end_timer = clock();
  std::cout << " Collapse Pyramid Cost " << (float)(end_timer - start_timer) / CLOCKS_PER_SEC * 1000 << " ms" <<  std::endl;
#endif
    return 0;
}


//==============

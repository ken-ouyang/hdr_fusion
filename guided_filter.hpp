#ifndef _GUIDED_FILTER_HPP_
#define _GUIDED_FILTER_HPP_


#include <opencv2/core/core.hpp>

using namespace cv;

void guidedfilter(const Mat &guided, const Mat &src, Mat &dst, int r, float eps);
void fastguidedfilter(const Mat &guided, const Mat &src, Mat &dst, int r, float eps, int s);
void fastguidedfilterColor(const Mat &guided, const Mat &src, Mat &dst, int r, float eps, int s);

#endif

#include "guided_filter.hpp"
#include <opencv2/opencv.hpp>

// Guided Filter for gray-scale image (Mat type: CV_32F)
void guidedfilter(const Mat &guided, const Mat &src, Mat &dst, int r, float eps)
{
    int H = guided.rows, W = guided.cols;

    Size ksize(r*2+1, r*2+1);
    Mat N, mean_I, mean_p, cov_Ip, var_I;
    const Mat &I = guided, &p = src;
    Mat q;
    Mat Ip = I.mul(p), II = I.mul(I);

    boxFilter(Mat::ones(H, W, guided.type()), N, -1, ksize, Point(-1, -1), false);
    boxFilter(I,  mean_I,  -1, ksize, Point(-1, -1), false); mean_I /= N;
    boxFilter(p,  mean_p,  -1, ksize, Point(-1, -1), false); mean_p /= N;
    boxFilter(Ip, cov_Ip, -1, ksize, Point(-1, -1), false); cov_Ip /= N;
    cov_Ip = cov_Ip - mean_I.mul(mean_p);

    boxFilter(II, var_I, -1, ksize, Point(-1, -1), false); var_I /= N;
    var_I = var_I - mean_I.mul(mean_I);

    Mat a, b;
    divide(cov_Ip, var_I + eps, a);     cov_Ip.release(); var_I.release();
    b = mean_p - a.mul(mean_I);         mean_I.release(); mean_p.release();

    //Mat mean_a, mean_b;
    boxFilter(a,  a,  -1, ksize, Point(-1, -1), false); a /= N;
    boxFilter(b,  b,  -1, ksize, Point(-1, -1), false); b /= N;

    q = a.mul(I) + b;
    q.copyTo(dst);
}

// FAST Guided Filter for gray-scale image (Mat type: CV_32F)
void fastguidedfilter(const Mat &guided, const Mat &src, Mat &dst, int r, float eps, int s)
{
    if(r < s)
		r = s;

    int H = guided.rows, W = guided.cols;

    int H_sub = H/s;
    int W_sub = W/s;
    int r_sub = r/s;
    Size size_sub(W_sub,H_sub);
    Size size(W,H);

    Mat guided_sub, src_sub;

    resize(guided,guided_sub,size_sub);
    resize(src,src_sub,size_sub);

    Size ksize(r_sub*2+1, r_sub*2+1);
    Mat N, mean_I, mean_p, cov_Ip, var_I;
    const Mat &I = guided_sub, &p = src_sub;
    Mat Ip = I.mul(p), II = I.mul(I);

    boxFilter(Mat::ones(H_sub, W_sub, guided_sub.type()), N, -1, ksize, Point(-1, -1), false);
    boxFilter(I,  mean_I,  -1, ksize, Point(-1, -1), false); mean_I /= N;

    boxFilter(p,  mean_p,  -1, ksize, Point(-1, -1), false); mean_p /= N;
    boxFilter(Ip, cov_Ip, -1, ksize, Point(-1, -1), false); cov_Ip /= N;  Ip.release();
    cov_Ip -= mean_I.mul(mean_p);

    boxFilter(II, var_I, -1, ksize, Point(-1, -1), false); var_I /= N;  II.release();
    var_I -= mean_I.mul(mean_I);

    Mat a_sub, b_sub;
    divide(cov_Ip, var_I + eps, a_sub);     cov_Ip.release(); var_I.release();
    b_sub = mean_p - a_sub.mul(mean_I);     mean_I.release(); mean_p.release();

    //Mat mean_a, mean_b;
    boxFilter(a_sub,  a_sub,  -1, ksize, Point(-1, -1), false); a_sub /= N;
    boxFilter(b_sub,  b_sub,  -1, ksize, Point(-1, -1), false); b_sub /= N;
    N.release();

    guided_sub.release();
    src_sub.release();

    Mat a, b;

    resize(a_sub,a,size);
    resize(b_sub,b,size);

    a_sub.release();
    b_sub.release();

    dst = a.mul(guided) + b;
	a.release();
	b.release();
}

// FAST Guided Filter for gray-scale image (Mat type: CV_32F), Has Bug, will fix it later. Temporally solution: run 3 channels independently.
// Created and Modified By yuwing@sensetime.com 7 Dec 2015
void fastguidedfilterColor(const Mat &guided, const Mat &src, Mat &dst, int r, float eps, int s)
{
    Mat guided_ch[3], src_ch[3], dst_ch[3];

    split(guided, guided_ch);
    split(src, src_ch);
    fastguidedfilter(guided_ch[0], src_ch[0], dst_ch[0], r, eps, s);
    fastguidedfilter(guided_ch[1], src_ch[1], dst_ch[1], r, eps, s);
    fastguidedfilter(guided_ch[2], src_ch[2], dst_ch[2], r, eps, s);

    merge(dst_ch, 3, dst);
/*
  if(r < s)
  r = s;

  int H = guided.rows, W = guided.cols;

  int H_sub = H/s;
  int W_sub = W/s;
  int r_sub = r/s;
  Size size_sub(W_sub,H_sub);
  Size size(W,H);

  Mat guided_ch[3], src_ch[3];
  Mat guided_sub_ch[3], src_sub_ch[3];

  split(guided, guided_ch);
  split(src, src_ch);

  resize(guided_ch[0],guided_sub_ch[0],size_sub);
  resize(guided_ch[1],guided_sub_ch[1],size_sub);
  resize(guided_ch[2],guided_sub_ch[2],size_sub);
  resize(src_ch[0],src_sub_ch[0],size_sub);
  resize(src_ch[1],src_sub_ch[1],size_sub);
  resize(src_ch[2],src_sub_ch[2],size_sub);

  Size ksize(r_sub*2+1, r_sub*2+1);
  Mat N, mean_I_r, mean_I_g, mean_I_b, mean_p_r, mean_p_g, mean_p_b, cov_Ip_r, cov_Ip_g, cov_Ip_b, var_I_rr, var_I_rg, var_I_rb, var_I_gg, var_I_gb, var_I_bb;

  Mat *I_ch[3], *p_ch[3];
  I_ch[0] = &guided_sub_ch[0];
  I_ch[1] = &guided_sub_ch[1];
  I_ch[2] = &guided_sub_ch[2];
  p_ch[0] = &src_sub_ch[0];
  p_ch[1] = &src_sub_ch[1];
  p_ch[2] = &src_sub_ch[2];

  Mat q;
  Mat Ip[3] = {I_ch[0]->mul(*p_ch[0]), I_ch[1]->mul(*p_ch[1]), I_ch[2]->mul(*p_ch[2])};

  boxFilter(Mat::ones(H_sub, W_sub, guided_sub_ch[0].type()), N, -1, ksize, Point(-1, -1), false);
  boxFilter(*I_ch[0],  mean_I_r,  -1, ksize, Point(-1, -1), false); mean_I_r /= N;
  boxFilter(*I_ch[1],  mean_I_g,  -1, ksize, Point(-1, -1), false); mean_I_g /= N;
  boxFilter(*I_ch[2],  mean_I_b,  -1, ksize, Point(-1, -1), false); mean_I_b /= N;

  boxFilter(*p_ch[0],  mean_p_r,  -1, ksize, Point(-1, -1), false); mean_p_r /= N;
  boxFilter(*p_ch[1],  mean_p_g,  -1, ksize, Point(-1, -1), false); mean_p_g /= N;
  boxFilter(*p_ch[2],  mean_p_b,  -1, ksize, Point(-1, -1), false); mean_p_b /= N;

  boxFilter(Ip[0], cov_Ip_r, -1, ksize, Point(-1, -1), false); cov_Ip_r /= N;
  boxFilter(Ip[1], cov_Ip_g, -1, ksize, Point(-1, -1), false); cov_Ip_g /= N;
  boxFilter(Ip[2], cov_Ip_b, -1, ksize, Point(-1, -1), false); cov_Ip_b /= N;
  cov_Ip_r = cov_Ip_r - mean_I_r.mul(mean_p_r);
  cov_Ip_g = cov_Ip_g - mean_I_g.mul(mean_p_g);
  cov_Ip_b = cov_Ip_b - mean_I_b.mul(mean_p_b);

  boxFilter(I_ch[0]->mul(*I_ch[0]), var_I_rr, -1, ksize, Point(-1,-1), false); var_I_rr /= N; var_I_rr = var_I_rr - mean_I_r.mul(mean_I_r) + eps;
  boxFilter(I_ch[0]->mul(*I_ch[1]), var_I_rg, -1, ksize, Point(-1,-1), false); var_I_rg /= N; var_I_rg = var_I_rg - mean_I_r.mul(mean_I_g);
  boxFilter(I_ch[0]->mul(*I_ch[2]), var_I_rb, -1, ksize, Point(-1,-1), false); var_I_rb /= N; var_I_rb = var_I_rb - mean_I_r.mul(mean_I_b);
  boxFilter(I_ch[1]->mul(*I_ch[1]), var_I_gg, -1, ksize, Point(-1,-1), false); var_I_gg /= N; var_I_gg = var_I_gg - mean_I_g.mul(mean_I_g) + eps;
  boxFilter(I_ch[1]->mul(*I_ch[2]), var_I_gb, -1, ksize, Point(-1,-1), false); var_I_gb /= N; var_I_gb = var_I_gb - mean_I_g.mul(mean_I_b);
  boxFilter(I_ch[2]->mul(*I_ch[2]), var_I_bb, -1, ksize, Point(-1,-1), false); var_I_bb /= N; var_I_bb = var_I_bb - mean_I_b.mul(mean_I_b) + eps;

  Mat invrr, invrg, invrb, invgg, invgb, invbb;
  // Inverse of Sigma + eps * I
  invrr = var_I_gg.mul(var_I_bb) - var_I_gb.mul(var_I_gb);
  invrg = var_I_gb.mul(var_I_rb) - var_I_rg.mul(var_I_bb);
  invrb = var_I_rg.mul(var_I_gb) - var_I_gg.mul(var_I_rb);
  invgg = var_I_rr.mul(var_I_bb) - var_I_rb.mul(var_I_rb);
  invgb = var_I_rb.mul(var_I_rg) - var_I_rr.mul(var_I_gb);
  invbb = var_I_rr.mul(var_I_gg) - var_I_rg.mul(var_I_rg);

  Mat covDet = invrr.mul(var_I_rr) + invrg.mul(var_I_rg) + invrb.mul(var_I_rb);

  invrr /= covDet;
  invrg /= covDet;
  invrb /= covDet;
  invgg /= covDet;
  invgb /= covDet;
  invbb /= covDet;

  Mat a_sub_ch[3], b_sub_ch[3];
  a_sub_ch[0] = invrr.mul(cov_Ip_r) + invrg.mul(cov_Ip_g) + invrb.mul(cov_Ip_b);
  a_sub_ch[1] = invrg.mul(cov_Ip_r) + invgg.mul(cov_Ip_g) + invgb.mul(cov_Ip_b);
  a_sub_ch[2] = invrb.mul(cov_Ip_r) + invgb.mul(cov_Ip_g) + invbb.mul(cov_Ip_b);

  b_sub_ch[0] = mean_p_r - a_sub_ch[0].mul(mean_I_r) - a_sub_ch[1].mul(mean_I_g) - a_sub_ch[2].mul(mean_I_b);
  b_sub_ch[1] = mean_p_g - a_sub_ch[0].mul(mean_I_r) - a_sub_ch[1].mul(mean_I_g) - a_sub_ch[2].mul(mean_I_b);
  b_sub_ch[2] = mean_p_b - a_sub_ch[0].mul(mean_I_r) - a_sub_ch[1].mul(mean_I_g) - a_sub_ch[2].mul(mean_I_b);

  //Mat mean_a, mean_b;
  boxFilter(a_sub_ch[0],  a_sub_ch[0],  -1, ksize, Point(-1, -1), false); a_sub_ch[0] /= N;
  boxFilter(a_sub_ch[1],  a_sub_ch[1],  -1, ksize, Point(-1, -1), false); a_sub_ch[1] /= N;
  boxFilter(a_sub_ch[2],  a_sub_ch[2],  -1, ksize, Point(-1, -1), false); a_sub_ch[2] /= N;
  boxFilter(b_sub_ch[0],  b_sub_ch[0],  -1, ksize, Point(-1, -1), false); b_sub_ch[0] /= N;
  boxFilter(b_sub_ch[1],  b_sub_ch[1],  -1, ksize, Point(-1, -1), false); b_sub_ch[1] /= N;
  boxFilter(b_sub_ch[2],  b_sub_ch[2],  -1, ksize, Point(-1, -1), false); b_sub_ch[2] /= N;

  Mat a_ch[3], b_ch[3], q_ch[3];

  resize(a_sub_ch[0],a_ch[0],size);
  resize(a_sub_ch[1],a_ch[1],size);
  resize(a_sub_ch[2],a_ch[2],size);
  resize(b_sub_ch[0],b_ch[0],size);
  resize(b_sub_ch[1],b_ch[1],size);
  resize(b_sub_ch[2],b_ch[2],size);

  q_ch[0] = a_ch[0].mul(guided_ch[0]) + b_ch[0];
  q_ch[1] = a_ch[1].mul(guided_ch[1]) + b_ch[0];
  q_ch[2] = a_ch[2].mul(guided_ch[2]) + b_ch[0];

  merge(q_ch, 3, dst);*/
}

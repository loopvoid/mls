#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

#ifdef DEBUG
#pragma comment(lib, "opencv_world341d.lib")
#else
#pragma comment(lib, "opencv_world341.lib")
#endif // DEBUG

Mat Gauss(double beta, Mat& r) {
	double b2 = beta * beta;
	Mat r2;
	multiply(r, r, r2);
	double eb2 = exp(-b2);

	Mat w(r.size(), CV_64F);

	for (int i = 0; i < w.rows; i++){
		for (int j = 0; j < w.cols; j++){
			double r_item = r.at<double>(i, j);
			double r2_item = r2.at<double>(i, j);
			if (r_item<=1.0) {
				w.at<double>(i, j) = (exp(-b2 * r2_item) - eb2) / (1.0 - eb2);
			}
			else{
				w.at<double>(i, j) = 0.0;
			}	
		}
	}
	return w;
}

/*
type : 权函数类型。Type = ("GAUSS", "CUBIC", "SPLI3", "SPLI5", "POWER", "CRBF1", "CRBF2", "CRBF3", "CRBF4", "CRBF5", "CRBF6")
para : 权函数参数
x    : 待评估其MLS形函数的点的x坐标
y    : 待评估其MLS形函数的点的y坐标
xI   : 用于构造MLS拟合函数的节点x坐标
yI   : 用于构造MLS拟合函数的节点y坐标
radius: 节点的支持域(影响域/影响半径)
*/
Mat rectangleWeight(string type, double para, int x, int y, Mat& xI, Mat& yI, double radius) {
	Mat rx = Mat::zeros(xI.size(), xI.type());
	Mat ry = Mat::zeros(yI.size(), yI.type());

	for (int i = 0; i < xI.rows; i++)
		for (int j = 0; j < xI.cols; j++)
			rx.at<double>(i, j) = abs(x - xI.at<double>(i, j)) / radius;

	for (int i = 0; i < yI.rows; i++)
		for (int j = 0; j < yI.cols; j++)
			ry.at<double>(i, j) = abs(y - yI.at<double>(i, j)) / radius;
	
	Mat wx, wy;
	if (type=="GAUSS")
	{
		wx = Gauss(para, rx);
		wy = Gauss(para, ry);
	}
	Mat w;
	multiply(wx, wy, w);
	return w;
}

Mat mls2Dshape(int m, int nnodes, Mat& xI, Mat& yI, int img_rows, int img_cols,
	double radius, string type, double para) {

	Mat DmI, wII;
	Mat wI = Mat::zeros(nnodes, nnodes, CV_64F);
	Mat xII = Mat::zeros(1, nnodes, CV_64F);
	Mat yII = Mat::zeros(1, nnodes, CV_64F);
	Mat B, pxy, B_T;
	int npoints = img_rows * img_cols;
	Mat PHI = Mat::zeros(npoints, nnodes, CV_64F);

	for (size_t ii = 0; ii < img_cols; ii++)
	{
		for (size_t jj = 0; jj < img_rows; jj++)
		{
			wII = rectangleWeight(type, para, ii, jj, xI, yI, radius);

			//wI=diag(wII(:));
			for (size_t i = 0; i < wII.cols; i++) {
				for (size_t j = 0; j < wII.rows; j++) {
					int diag_index = i * wII.rows + j;
					wI.at<double>(diag_index, diag_index) = wII.at<double>(j, i);
				}
			}

			//xII=reshape(xI,size(xII));
			for (size_t i = 0; i < xI.cols; i++)
				for (size_t j = 0; j < xI.rows; j++)
					xII.at<double>(0, i * xI.rows + j) = xI.at<double>(j, i);

			//yII=reshape(yI,size(yII));
			for (size_t i = 0; i < yI.cols; i++)
				for (size_t j = 0; j < yI.rows; j++)
					yII.at<double>(0, i * yI.rows + j) = yI.at<double>(j, i);

			if (m == 6) {
				Mat xII2, xIIyII, yII2;
				multiply(xII, xII, xII2);
				multiply(xII, yII, xIIyII);
				multiply(yII, yII, yII2);

				//B = [ones(1, nnodes); xII;yII; xII.*xII;xII.*yII;yII.*yII]; 
				Mat B_temp(m, nnodes, CV_64F);
				for (size_t i = 0; i < nnodes; i++) {
					B_temp.at<double>(0, i) = 1;
					B_temp.at<double>(1, i) = xII.at<double>(0, i);
					B_temp.at<double>(2, i) = yII.at<double>(0, i);
					B_temp.at<double>(3, i) = xII2.at<double>(0, i);
					B_temp.at<double>(4, i) = xIIyII.at<double>(0, i);
					B_temp.at<double>(5, i) = yII2.at<double>(0, i);
				}
				B = B_temp;

				//pxy = [1; x(j); y(j);x(j)*x(j);x(j)*y(j);y(j)*y(j)];
				pxy = Mat(6, 1, CV_64F);
				pxy.at<double>(0, 0) = 1;
				pxy.at<double>(1, 0) = ii;
				pxy.at<double>(2, 0) = jj;
				pxy.at<double>(3, 0) = ii * ii; 
				pxy.at<double>(4, 0) = ii * jj;
				pxy.at<double>(5, 0) = jj * jj;
				//pxy = Mat(6, 1, CV_64F);
				//pxy.at<float>(0, 0) = 1;
				//pxy.at<float>(1, 0) = x_index;
				//pxy.at<float>(2, 0) = y_index;
				//pxy.at<float>(3, 0) = x_index * y_index;
				//pxy.at<float>(4, 0) = x_index * x_index;
				//pxy.at<float>(5, 0) = y_index * y_index;
			}

			cv::transpose(B, B_T);
			Mat BwI_temp = B * wI;
			Mat A = BwI_temp * B_T;
			//Mat A = B * wI * B_T;
			
			double temp_A = A.at<double>(A.rows - 1, A.cols - 1);

			double radius_plus = radius;

			Mat A_inv = Mat::zeros(A.size(),A.type());
			double RAcond = 0;
			RAcond = invert(A, A_inv, DECOMP_SVD);
						
			//while (RAcond <= 9.999999e-016)
			int counter = 0;
			//while (RAcond <= 9.999999e-016)
			while (false)
			{
				counter++;
				if (counter > 2) break;
				radius_plus *= 2;
				wII = rectangleWeight(type, para, ii, jj, xI, yI, radius_plus);
				//wII = rectangleWeight(type, para, x_index, y_index, xI, yI, radius_plus);

				//wI=diag(wII(:));
				for (size_t i = 0; i < wII.cols; i++) {
					for (size_t j = 0; j < wII.rows; j++) {
						int diag_index = i * wII.rows + j;
						wI.at<double>(diag_index, diag_index) = wII.at<double>(j, i);
					}
				}

				A = B * wI * B_T;
				RAcond = invert(A, A_inv, DECOMP_SVD);
			}

			Mat pxy_T;
			transpose(pxy, pxy_T);
			Mat PHI_one_row = pxy_T * A_inv * B * wI;
			for (size_t i = 0; i < PHI_one_row.cols; i++) {
				double* p_PHI_DATA = PHI.ptr<double>(ii * img_rows + jj);
				double* p_PHI_one_row = PHI_one_row.ptr<double>(0);
				p_PHI_DATA[i] = p_PHI_one_row[i];
				//PHI.at<double>((ii-1)*img_rows+(jj-1), i) = PHI_one_row.at<double>(0, i);
				//PHI.at<double>(ii * img_rows + jj, i) = PHI_one_row.at<double>(0, i);
			}
				
		}
	}

	return PHI;
}
Mat mls2d(Mat& img, double radius=30, double para=3.0, string type="GAUSS") {
	//double radius = 30;
	//string type = "GAUSS";
	//double para = 3.0;

	int step = 10;
	int nodal_rows = (int)img.rows / step + 1;
	int nodal_cols = (int)img.cols / step + 1;
	int nnodes = nodal_cols * nodal_rows;
	int img_rows = img.rows;
	int img_cols = img.cols;
	int npoints = img_rows * img.cols;

	Mat xII(1, nodal_rows, CV_64F);
	Mat yII(1, nodal_cols, CV_64F);
	for (size_t i = 0; i < nodal_rows; i++)
		xII.at<double>(0, i) = i * 10 + 1;
	for (size_t i = 0; i < nodal_cols; i++)
		yII.at<double>(0, i) = i * 10 + 1;

	Mat xI(nodal_rows, nodal_cols, CV_64F);
	Mat yI(nodal_rows, nodal_cols, CV_64F);
	for (size_t i = 0; i < nodal_rows; i++)
		for (size_t j = 0; j < nodal_cols; j++)
			xI.at<double>(i, j) = j * 10 + 1;
	for (size_t i = 0; i < nodal_cols; i++)
		for (size_t j = 0; j < nodal_rows; j++)
			yI.at<double>(j, i) = j * 10 + 1;

	Mat x(img_rows, img_cols, CV_64F);
	Mat y(img_rows, img_cols, CV_64F);
	for (size_t i = 0; i < img_cols; i++)
		for (size_t j = 0; j < img_rows; j++)
			x.at<double>(j, i) = i;
	//for (size_t i = 0; i < img_rows; i++) x.at<double>(i, img_cols-1) = img_cols - 1;


	for (size_t i = 0; i < img_rows; i++)
		for (size_t j = 0; j < img_cols; j++)
			y.at<double>(i, j) = i;
	//for (size_t i = 0; i < img_cols; i++) y.at<double>(img_rows-1, i) = img_rows - 1;


	//Mat dmI(nodal_rows, nodal_cols, CV_64F);
	//for (size_t i = 0; i < nodal_rows; i++)
	//	for (size_t j = 0; j < nodal_cols; j++)
	//		dmI.at<float>(i, j) = radius;
		

	//Mat PHI = mls2Dshape(6, nnodes, xI, yI, img.rows, img.cols, x, y, radius, type, para);
	Mat PHI = mls2Dshape(6, nnodes, xI, yI, img.rows, img.cols, radius, type, para);


	Mat ZII(nodal_rows, nodal_cols, CV_64F);
	for (size_t i = 0; i < nodal_rows; i++)
		for (size_t j = 0; j < nodal_cols; j++)
			ZII.at<double>(i, j) = img.at<uchar>(xII.at<double>(0, i), yII.at<double>(0, j));

	Mat Znodes(1, nnodes, CV_64F);
	for (size_t i = 0; i < nodal_cols; i++)
		for (size_t j = 0; j < nodal_rows; j++)
			Znodes.at<double>(0, i * nodal_rows + j) = ZII.at<double>(j, i);
	
	Mat zh(npoints, 1, CV_64F);
	Mat Znodes_T;
	cv::transpose(Znodes, Znodes_T);
	zh = PHI * Znodes_T;

	Mat out(img_rows, img_cols, CV_64F);
	for (size_t i = 0; i < img_cols; i++)
		for (size_t j = 0; j < img_rows; j++)
			out.at<double>(j, i) = zh.at<double>(i * img_rows + j, 0);

	return out;
}

Mat removeAllzeroRowCol(Mat& in) {
	Mat src = in;
	Point left_top;
	left_top.x = 0;
	left_top.y = 0;
	int counter = 0;
	int thread = 5;
	for (size_t i = 0; i < src.rows; i++)
	{
		uchar* pData = src.ptr<uchar>(i);
		for (size_t j = 0; j < src.cols; j++)
		{
			if (pData[j] <=thread) counter++;
		}
		if (counter == src.cols) {
			left_top.y++;
			counter = 0;
		}
		else
			break;
		counter = 0;
	}

	counter = 0;
	for (size_t i = 0; i < src.cols; i++)
	{
		uchar* pData = src.ptr<uchar>(i);
		for (size_t j = 0; j < src.rows; j++)
		{
			if (src.at<uchar>(j, i) <= thread) counter++;
		}
		if (counter == src.rows) {
			left_top.x++;
			counter = 0;
		}
		else
			break;
		counter = 0;
	}

	Mat _roi = src(Rect(left_top, Point(src.cols, src.rows)));
	return _roi;
}

int main()
{
	double radius = 30;
	double para = 3.0;
	string type = "GAUSS";
	
	//Mat in = imread("img_8uc1.jpg", 0);
	//Mat roi = removeAllzeroRowCol(in);
	//imwrite("roi.jpg", roi);

	//Mat in = imread("dispmat.jpg", 0);
	//Mat in = imread("roi.jpg", 0);
	Mat in = imread("22result_ch1.jpg", 0);
	//resize(in, in, Size(in.cols/2 , in.rows/2 ));

	clock_t start = clock();
	Mat out = mls2d(in, radius, para, type);
	cout << "Time: " << clock() - start << endl;

	//resize(out, out, Size(in.cols*2 , in.rows*2 ));

	imshow("out", out);
	cvWaitKey(0);
	imwrite("aaaa.jpg", out);
	
	return 0;
}
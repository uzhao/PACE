// must enable EIGEN2_SUPPORT otherwise we can't return MatrixXd in some case
// should be removed after they fix this problem
// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2012-May/003792.html
// commoned variable name in parentheses is the old name in matlab version

#define EIGEN2_SUPPORT
#define GAUSSIAN_MULT 0.3989422804014326779399460599343818684758586311649346
#define LAMBDA_THRESHOLD 0.000001

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Generic Optimize Function for Class
// ax:  lower bound
// bx:  upper bound
// obj: pointer to the object
// f:   pointer to the method(function) which need optimize
//      input is a double(e.g.:bandwidth), output is a double(e.g.:gcv)
template<class T>
double gopt(double ax, double bx, T* obj, double(T::*f) (double)) {
	const double c = (3. - sqrt(5.)) * .5;
	double a, b, d, e, p, q, r, u, v, w, x;
	double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
	eps = DBL_EPSILON;
	tol1 = eps + 1.;
	eps = sqrt(eps);
	a = ax;
	b = bx;
	v = a + c * (b - a);
	w = v;
	x = v;
	d = 0.;
	e = 0.;
	fx = (obj->*f)(x);
	fv = fx;
	fw = fx;
	tol3 = (1e-9) / 3.;
	for (;;) {
		xm = (a + b) * .5;
		tol1 = eps * fabs(x) + tol3;
		t2 = tol1 * 2.;
		if (fabs(x - xm) <= t2 - (b - a) * .5) break;
		p = 0.;
		q = 0.;
		r = 0.;
		if (fabs(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = (q - r) * 2.;
			if (q > 0.) p = -p; else q = -q;
			r = e;
			e = d;
		}
		if (fabs(p) >= fabs(q * .5 * r) ||
			p <= q * (a - x) || p >= q * (b - x)) {
			if (x < xm) e = b - x; else e = a - x;
			d = c * e;
		}
		else {
			d = p / q;
			u = x + d;
			if (u - a < t2 || b - u < t2) {
				d = tol1;
				if (x >= xm) d = -d;
			}
		}
		if (fabs(d) >= tol1)
			u = x + d;
		else if (d > 0.)
			u = x + tol1;
		else
			u = x - tol1;
		fu = (obj->*f)(u);
		if (fu <= fx) {
			if (u < x) b = x; else a = x;
			v = w;    w = x;   x = u;
			fv = fw; fw = fx; fx = fu;
		}
		else {
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x) {
				v = w; fv = fw;
				w = u; fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u; fv = fu;
			}
		}
	}
	return x;
}

// Trapezoidal Numerical Integration
// Integration for f(x) = y
double trapz(VectorXd &x, VectorXd &y) {
	int n = x.size() - 1;
	return (x.segment(1, n) - x.segment(0, n)).cwiseProduct(y.segment(1, n) + y.segment(0, n)).sum() / 2;
}

// for mean, eigen func
// also support other Xlwls
class lwls {
public:
	VectorXd     grid;        // output grid (xou in lwls.m)
	VectorXd     output;      // output value (mu in lwls.m)
	VectorXd     new_weight;  // new weight for 2d smoothing, may not useful depends on smoothing method
	VectorXd     x;           // input grid (xin in lwls.m)
	VectorXd     y;           // input grid (yin in lwls.m)
	VectorXd     w;           // input grid (win in lwls.m)
	VectorXd     mean;        // output mean, for gcv
	VectorXd     hat;         // diag of hat matrix, for gcv
	int          degree;      // degree of polynomial (npoly in lwls.m)
	int          drv;         // order of derivative (nder in lwls.m)
	double       h;           // mean of diag of hat matrix, for gcv
	double       bw;          // bandwidth (bw in lwls.m)
	double       minbw;       // minimum bandwidth
	double       maxbw;       // maximum bandwidth
	double       gcvv;        // GCV value

	lwls(VectorXd x_, VectorXd y_, VectorXd w_, VectorXd grid_, int degree_, int drv_) : x(x_), y(y_), w(w_), grid(grid_), degree(degree_), drv(drv_) {
		output = VectorXd::Zero(grid.size());
		mean = VectorXd::Zero(grid.size());
		new_weight = VectorXd::Zero(grid.size());
		hat = VectorXd::Zero(grid.size());
		minbw = (x.tail(x.size() - 1) - x.head(x.size() - 1)).mean();
		maxbw = x[x.size() - 1] - x[0];
	};

	void optimize() {
		double best_gcv = gopt(minbw, maxbw, this, &lwls::update);
		bw = sqrt(best_gcv * minbw);  // using geometric mean for bandwidth
		update(bw, true);
	};

	// overload update with one agreement so we can apply gopt
	double update(double bw_) {
		return update(bw_, false);
	};

	// refer to 3rd ppt from Cong
	// if everything is true, then also output the new weight for smoothing on second direction
	double update(double bw_, bool everything){
		VectorXd true_grid(grid);
		if (!everything) {
			grid = x;
		}
		bw = bw_;
		VectorXd d = VectorXd::Zero(x.size());
		VectorXd d_bw = VectorXd::Zero(x.size());
		VectorXd k_d_bw = VectorXd::Zero(x.size());
		VectorXd w_k_d_bw = VectorXd::Zero(x.size());

		MatrixXd local_x;
		VectorXd local_respond;
		VectorXd design_base;
		VectorXd respond_base;
		VectorXd local_ans;
		VectorXd tempweight;

		for (int i = 0; i != grid.size(); i++) {
			d = -x;
			d.array() += grid[i];
			// Guassian kernel only
			d_bw = d / bw;
			for (int j = 0; j != x.size(); j++) {
				k_d_bw[j] = GAUSSIAN_MULT * exp(-0.5 * d_bw[j] * d_bw[j]);
			}
			w_k_d_bw = k_d_bw.cwiseProduct(w);

			if (everything) {
				tempweight = k_d_bw.array().square();
				new_weight[i] = tempweight.cwiseProduct(w).sum();
				if (new_weight[i] < DBL_EPSILON) {
					new_weight[i] = 0;
				}
				else {
					new_weight[i] = 1 / new_weight[i];
				}
			}

			local_x = MatrixXd::Zero(degree + 1, degree + 1);
			local_respond = VectorXd::Zero(degree + 1);
			design_base = w_k_d_bw;
			respond_base = w_k_d_bw.cwiseProduct(y);

			// fill first row
			for (int j = 0; j != degree + 1; j++) {
				local_x(0, j) = design_base.sum();
				local_respond[j] = respond_base.sum();
				design_base = design_base.cwiseProduct(d);
				respond_base = respond_base.cwiseProduct(d);
			}
			// fill last col, no more respond need to fill
			for (int j = 1; j != degree + 1; j++) {
				local_x(degree, j) = design_base.sum();
				design_base = design_base.cwiseProduct(d);
			}
			// fill the rest
			for (int row = 1; row != degree + 1; row++) {
				for (int col = 0; col != degree; col++) {
					local_x(row, col) = local_x(row - 1, col + 1);
				}
			}
			// FIXME: THIS HAT MATRIX IS ONLY FOR DEGREE 1
			hat[i] = w[i] * local_x(1, 1) / (local_x(0, 0)*local_x(1, 1) - local_x(0, 1)*local_x(1, 0));

			local_ans = local_x.jacobiSvd(ComputeThinU | ComputeThinV).solve(local_respond);
			mean[i] = local_ans[0];
			output[i] = local_ans[drv];
		}
		h = hat.mean(); // GCV method, remove this line if use CV
		gcvv = (y - mean).array().square().sum() / (1 - 2 * h + h * h);
		if (!everything) {
			grid = true_grid;
		}
		return gcvv;
	};

	~lwls() {
	};
};

// support other Xlwls
// smooth along column
class mlwls {
public:
	double             bw;
	double             minbw;
	double             maxbw;
	int                degree;
	int                drv;
	double             gcvv;
	int                kernel;
	std::vector<lwls*> lwlss;
	MatrixXd           smtmat;
	MatrixXd           new_weight;

	mlwls(VectorXd x, MatrixXd mat, MatrixXd weight, int degree_, int drv_) :
		degree(degree_), drv(drv_) {
		smtmat = MatrixXd::Zero(mat.rows(), mat.cols());
		new_weight = MatrixXd::Zero(mat.rows(), mat.cols());
		int n = mat.cols();
		lwlss = std::vector<lwls*>(n);
		minbw = (x.tail(x.size() - degree) - x.head(x.size() - degree)).maxCoeff();
		maxbw = (x.tail(1) - x.head(1))[0] / 2;
		for (int i = 0; i != n; i++) {
			lwlss[i] = new lwls(x, mat.col(i), weight.col(i), x, degree, drv);
		}
		optimize();
	};

	void optimize() {
		double best_gcv;
		best_gcv = gopt(minbw, maxbw, this, &mlwls::update);
		// bw = best_gcv;  // GCV
        // bw = minbw;     // bin grid
        bw = sqrt(best_gcv * minbw);  // geo mean
		update(bw, true);
	};

	double update(double bw_) {
		return update(bw_, false);
	};

	double update(double bw_, bool everything) {
		bw = bw_;
		gcvv = 0;
		for (int i = 0; i != lwlss.size(); i++) {
			gcvv += lwlss[i]->update(bw, everything);
			if (everything) {
				smtmat.col(i) = lwlss[i]->output;
				new_weight.col(i) = lwlss[i]->new_weight;
			}
		}
		return gcvv;
	};

	~mlwls() {
		for (std::vector<lwls*>::iterator it = lwlss.begin(); it != lwlss.end(); it++) {
			delete *it;
		}
	};
};

// for cov mat
class covlwls {
public:
	MatrixXd midmat;
	MatrixXd midweight;
	MatrixXd smtmat;

	covlwls(VectorXd x, MatrixXd mat, MatrixXd weight) {
		mlwls middle(x, mat, weight, 1, 0);
		midmat = middle.smtmat.transpose();
		midweight = middle.new_weight.transpose();
		// midweight = weight.transpose();

		mlwls final(x, midmat, midweight, 1, 0);
		smtmat = (final.smtmat + final.smtmat.transpose()) / 2;
	};

	~covlwls() {
	};
};

// diagonal elements of covariance matrix
// FIXME
//////////////////////////////////////////////////////////
class qlwls {
public:
	VectorXd smtdiag;

	qlwls(MatrixXd mat, MatrixXd weight, int boundary) {
		mat = mat.rowwise().reverse();
		weight = weight.rowwise().reverse();
	};

	~qlwls() {
	};
};

// for partial cov mat
// can not be used until fix hat matrix in lwls and
// figure out is bandwidth and so on
// just copy from covlwls with new degree and drv
// class pcovlwls {
// public:
// };

// class for individual subject
class Subject {
public:
	VectorXd x;
	VectorXd y;
	// VectorXd binx;
	VectorXd biny;

	bool     valid;  // if it has no observation in cutted domain, then valid is false
	VectorXd validx;
	VectorXd validy;

	VectorXd *invmean; // pointer to its fitted mean
	MatrixXd *invcov;  // pointer to its fitted cov matrix

	VectorXd pcs; // pc scores

	// irregular only
	VectorXd count;

	Subject(VectorXd &x_, VectorXd &y_) : x(x_), y(y_){
	};

	~Subject(){
		if (invmean != NULL) {
			delete invmean;
			invmean = NULL;
			delete invcov;
			invcov = NULL;
		}
	};
};

class FPCAreg {
public:
	int      n;            // number of subjects
	int      b;            // number of bins
	int      c;            // number of components

	std::vector<Subject*> subjects;

	VectorXd grid;
	VectorXd cuttedgrid;

	VectorXd cross_sectional_mean;
	VectorXd mean_weight;
	VectorXd mean;
	VectorXd cuttedmean;

	MatrixXd rawcov;
	MatrixXd cov_weight;
	VectorXd diag_weight;
	MatrixXd smtcov;
	MatrixXd cuttedsmtcov;
	VectorXd lambda;
	MatrixXd eigenfunc;
	MatrixXd cuttedfitcov;

	MatrixXd pcs;

	double   cutp;
	int      boundary;
	int      cuttedlength;
	double   minx;
	double   maxx;
	double   cutminx;
	double   cutmaxx;
	double   gap;
	VectorXd fve;
	double   fveth;
	int      posk;
	int      maxk;

	bool     error;
	bool     new_ridge;
	double   sigma;
	double   rho;

	void set_smt_cov() {
		covlwls covlwls_cov = covlwls(grid, rawcov, cov_weight);
		smtcov = covlwls_cov.smtmat;
		// FIXME: fix diagonal elements of covmat
		// qlwls qlwls_diag = qlwls(rawcov, cov_weight, boundary);
		// smtcov.diagonal() = qlwls_diag.smtdiag;
		cuttedsmtcov = smtcov.block(boundary, boundary, cuttedlength, cuttedlength);
	};

	void set_eigen() {
		// eigen decomposition
		SelfAdjointEigenSolver<MatrixXd> eigendcp(cuttedsmtcov);
		lambda = eigendcp.eigenvalues().reverse();
		posk = 0;
		double lambdasum = (lambda.sum() + lambda.cwiseAbs().sum()) / 2;
		for (int i = 0; i != lambda.size(); i++) {
			if (lambda[i] / lambdasum > LAMBDA_THRESHOLD) {
				posk++;
			}
		}
		lambda = eigendcp.eigenvalues().reverse().head(posk) * gap;
		// choose number of components
		c = 2;
		fve = lambda / lambda.sum();
		for (int i = 1; i != fve.size(); i++) {
			fve[i] += fve[i - 1];
			if (fve[i] < fveth) {
				c = i + 2;
			}
		}
		c = c > maxk ? maxk : c;
		// set eigen value and eigen function
		lambda = eigendcp.eigenvalues().reverse().head(c) * gap;
		eigenfunc = eigendcp.eigenvectors().rowwise().reverse().leftCols(c) / sqrt(gap);
	};

	void set_fit_cov() {
		cuttedfitcov.setZero(cuttedlength, cuttedlength);
		for (int i = 0; i != c; i++) {
			cuttedfitcov += lambda[i] * eigenfunc.col(i) * eigenfunc.col(i).transpose();
		}
	};

	void set_sigma() {
		lwls covdiag = lwls(cuttedgrid, rawcov.diagonal().segment(boundary, cuttedlength), diag_weight.segment(boundary, cuttedlength), cuttedgrid, 1, 0);
		covdiag.optimize();
		VectorXd vare = covdiag.output - cuttedsmtcov.diagonal();
		sigma = trapz(cuttedgrid, vare) / (cuttedgrid[cuttedlength - 1] - cuttedgrid[0]);
		sigma = sigma > 0 ? sigma : 0;
	};

	void fpca() {
		set_mean();
		set_raw_cov();
		set_smt_cov();
		set_eigen();
		set_fit_cov();
		if (error) {
			set_sigma();
		}
		else {
			sigma = 0;
		}
		update_ridge();
		set_score();
	};

	FPCAreg(List x_, List y_, double cutp_, double fveth_, int maxk_, bool error_) :
		cutp(cutp_), fveth(fveth_), maxk(maxk_), error(error_) {
		grid = x_[0];                          // for regular case
		n = x_.size();                         // get number of subjects
		subjects = std::vector<Subject*>(n);   // initial subject points
		set_bin_no();                          // set number of bins
		set_subjects(x_, y_);                  // set subjects
		maxx *= 1 + DBL_EPSILON;               // to avoid overflow
		set_grids();
		cut_subjects();
		set_csm();
		fpca();
	};

	~FPCAreg() {
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			delete *it;
		}
	};

	void set_bin_no() {
		b = grid.size();
	};

	void set_subjects(List &x_, List &y_){
		maxx = grid.maxCoeff();
		minx = grid.minCoeff();
		for (int i = 0; i != n; i++) {
			VectorXd yt = y_[i];
			subjects[i] = new Subject(grid, yt);
		}
	};

	void set_grids() {
		gap = (maxx - minx) / (b - 1);
		boundary = ceil(b * cutp);
		cuttedlength = b - 2 * boundary;
		cuttedgrid = grid.segment(boundary, cuttedlength);
		cutminx = cuttedgrid.minCoeff();
		cutmaxx = cuttedgrid.maxCoeff();
	};

	void cut_subject(Subject &s) {
		// s.binx = s.x;
		s.biny = s.y;
		s.valid = true;
	};

	void cut_subjects() {
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			cut_subject(**it);
		}
	};

	void set_csm() { //set cross sectional mean
		cross_sectional_mean = VectorXd::Zero(b);
		mean_weight = VectorXd::Constant(b, n);
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			cross_sectional_mean += (*it)->y;
		}
		cross_sectional_mean /= n;
	};

	void set_mean() {
		mean = cross_sectional_mean;
	};

	void set_raw_cov() {
		rawcov.setZero(b, b);
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			rawcov += ((*it)->y - mean) * ((*it)->y - mean).transpose();
		}
		rawcov /= n;
		cov_weight.setConstant(b, b, n);
		diag_weight.setConstant(b, n);
		if (error) {
			cov_weight.diagonal().setZero();
		}
	};

	void update_ridge() {
		// TODO
		rho = sigma;
	};

	void set_score() {
		//
	};
};

class FPCAirreg {
public:
	int      n;            // number of subjects
	int      b;            // number of bins
	int      c;            // number of components

	std::vector<Subject*> subjects;

	VectorXd grid;
	VectorXd cuttedgrid;

	VectorXd cross_sectional_mean;
	VectorXd mean_weight;
	VectorXd mean;
	VectorXd cuttedmean;

	MatrixXd rawcov;
	MatrixXd cov_weight;
	VectorXd diag_weight;
	MatrixXd smtcov;
	MatrixXd cuttedsmtcov;
	VectorXd lambda;
	MatrixXd eigenfunc;
	MatrixXd cuttedfitcov;

	MatrixXd pcs;

	double   cutp;
	int      boundary;
	int      cuttedlength;
	double   minx;
	double   maxx;
	double   cutminx;
	double   cutmaxx;
	double   gap;
	VectorXd fve;
	double   fveth;
	int      posk;
	int      maxk;

	bool     error;
	bool     new_ridge;
	double   sigma;
	double   rho;

	void set_smt_cov() {
		covlwls covlwls_cov = covlwls(grid, rawcov, cov_weight);
		smtcov = covlwls_cov.smtmat;
		// FIXME: fix diagonal elements of covmat
		// qlwls qlwls_diag = qlwls(rawcov, cov_weight, boundary);
		// smtcov.diagonal() = qlwls_diag.smtdiag;
		cuttedsmtcov = smtcov.block(boundary, boundary, cuttedlength, cuttedlength);
	};

	void set_eigen() {
		// eigen decomposition
		SelfAdjointEigenSolver<MatrixXd> eigendcp(cuttedsmtcov);
		lambda = eigendcp.eigenvalues().reverse();
		posk = 0;
		double lambdasum = (lambda.sum() + lambda.cwiseAbs().sum()) / 2;
		for (int i = 0; i != lambda.size(); i++) {
			if (lambda[i] / lambdasum > LAMBDA_THRESHOLD) {
				posk++;
			}
		}
		lambda = eigendcp.eigenvalues().reverse().head(posk) * gap;
		// choose number of components
		c = 2;
		fve = lambda / lambda.sum();
		for (int i = 1; i != fve.size(); i++) {
			fve[i] += fve[i - 1];
			if (fve[i] < fveth) {
				c = i + 2;
			}
		}
		c = c > maxk ? maxk : c;
		// set eigen value and eigen function
		lambda = eigendcp.eigenvalues().reverse().head(c) * gap;
		eigenfunc = eigendcp.eigenvectors().rowwise().reverse().leftCols(c) / sqrt(gap);
	};

	void set_fit_cov() {
		cuttedfitcov.setZero(cuttedlength, cuttedlength);
		for (int i = 0; i != c; i++) {
			cuttedfitcov += lambda[i] * eigenfunc.col(i) * eigenfunc.col(i).transpose();
		}
	};

	void set_sigma() {
		lwls covdiag = lwls(cuttedgrid, rawcov.diagonal().segment(boundary, cuttedlength), diag_weight.segment(boundary, cuttedlength), cuttedgrid, 1, 0);
		covdiag.optimize();
		VectorXd vare = covdiag.output - cuttedsmtcov.diagonal();
		sigma = trapz(cuttedgrid, vare) / (cuttedgrid[cuttedlength - 1] - cuttedgrid[0]);
		sigma = sigma > 0 ? sigma : 0;
	};

	void fpca() {
		set_mean();
		set_raw_cov();
		set_smt_cov();
		set_eigen();
		set_fit_cov();
		if (error) {
			set_sigma();
		}
		else {
			sigma = 0;
		}
		update_ridge();
		set_score();
	};

	FPCAirreg(List x_, List y_, double cutp_, double fveth_, int maxk_, bool error_) :
		cutp(cutp_), fveth(fveth_), maxk(maxk_), error(error_) {
		// grid = x_[0];                          // for regular case
		n = x_.size();                         // get number of subjects
		subjects = std::vector<Subject*>(n);   // initial subject points
		set_bin_no();                          // set number of bins
		set_subjects(x_, y_);                  // set subjects
		maxx *= 1 + DBL_EPSILON;               // to avoid overflow
		set_grids();
		cut_subjects();
		set_csm();
		fpca();
	};

	~FPCAirreg() {
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			delete *it;
		}
	};

	void set_bin_no() {
		if (n >= 400) {
			b = 100;
		}
		else if (n <= 100) {
			b = 50;
		}
		else {
			b = 5 * sqrt(n);
		}
	};

	void set_subjects(List &x_, List &y_){
		for (int i = 0; i != n; i++) {
			VectorXd xt = x_[i];
			VectorXd yt = y_[i];
			if (i == 0) {
				maxx = xt.maxCoeff();
				minx = xt.minCoeff();
			}
			maxx = maxx > xt.maxCoeff() ? maxx : xt.maxCoeff();
			minx = minx < xt.minCoeff() ? minx : xt.minCoeff();
			subjects[i] = new Subject(xt, yt);
		}
	};

	void set_grids() {
		gap = (maxx - minx) / b;
		boundary = ceil(b * cutp);
		cuttedlength = b - 2 * boundary;
		grid = VectorXd::LinSpaced(b, minx + gap / 2, maxx - gap / 2);
		cuttedgrid = grid.segment(boundary, cuttedlength);
		cutminx = cuttedgrid.minCoeff();
		cutmaxx = cuttedgrid.maxCoeff();
	};

	void cut_subject(Subject &s) {
		s.count.setZero(b);
		s.biny.setZero(b);
		// s.binx = s.x;
		for (int i = 0; i != s.x.size(); i++) {
			int ii = (s.x[i] - minx) / gap;
			s.count[ii]++;
			s.biny[ii] += s.y[i];
		}
		for (int i = 0; i != b; i++) {
			if (s.count[i] != 0) {
				s.biny[i] /= s.count[i];
			}
		}
		s.valid = s.count.segment(boundary, cuttedlength).sum() > 0;
	};

	void cut_subjects() {
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			cut_subject(**it);
		}
	};

	void set_csm() { //set cross sectional mean
		cross_sectional_mean = VectorXd::Zero(b);
		mean_weight = VectorXd::Zero(b);
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			cross_sectional_mean += (*it)->biny;
			mean_weight = mean_weight + (*it)->count;
		}
		cross_sectional_mean = cross_sectional_mean.cwiseQuotient(mean_weight);
		for (int i = 0; i != b; i++) {
			if (mean_weight[i] == 0) {
				cross_sectional_mean[i] = 0;
			}
		}
	};

	void set_mean() {
		mean = cross_sectional_mean;
		// lwls lwls_mean = lwls(grid, cross_sectional_mean, mean_weight, grid, 1, 0);
		// lwls_mean.optimize();
		// mean = lwls_mean.output;
	};

	void set_raw_cov() {
		rawcov.setZero(b, b);
		cov_weight.setZero(b, b);
		double ei, ej;
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			for (int i = 0; i != b; i++) {
				if ((*it)->count[i] != 0) {
					for (int j = i; j != b; j++) {
						if ((*it)->count[j] != 0) {
							ei = (*it)->biny[i] - mean[i];
							ej = (*it)->biny[j] - mean[j];
							rawcov(i, j) += ei * ej;
							rawcov(j, i) += ei * ej;
							cov_weight(i, j)++;
							cov_weight(j, i)++;
						}
					}
				}
			}
		}
		rawcov.diagonal() /= 2;
		cov_weight.diagonal() /= 2;
		diag_weight = cov_weight.diagonal();
		for (int i = 0; i != b; i++) {
			for (int j = 0; j != b; j++) {
				if (cov_weight(i, j) >= 1) {
					rawcov(i, j) /= cov_weight(i, j);
				}
			}
		}
		if (error) {
			cov_weight.diagonal().setZero();
		}
	};

	void update_ridge() {
		// TODO
		rho = sigma;
	};

	void set_score() {
		//
	};
};


RCPP_MODULE(FPCA){
	using namespace Rcpp;

	class_<FPCAreg>("FPCAreg")
		.constructor<List, List, double, double, int, bool>()
		.field("n", &FPCAreg::n)
		.field("b", &FPCAreg::b)
		.field("c", &FPCAreg::c)
		.field("gap", &FPCAreg::gap)

		.field("grid", &FPCAreg::grid)
		.field("cuttedgrid", &FPCAreg::cuttedgrid)

		.field("cross_sectional_mean", &FPCAreg::cross_sectional_mean)
		.field("mean", &FPCAreg::mean)

		.field("rawcov", &FPCAreg::rawcov)
		.field("smtcov", &FPCAreg::smtcov)
		.field("cuttedsmtcov", &FPCAreg::cuttedsmtcov)
		.field("cuttedfitcov", &FPCAreg::cuttedfitcov)

		.field("lambda", &FPCAreg::lambda)
		.field("eigenfunc", &FPCAreg::eigenfunc)

		.field("sigma", &FPCAreg::sigma)

		;

	class_<FPCAirreg>("FPCAirreg")
		.constructor<List, List, double, double, int, bool>()
		.field("n", &FPCAirreg::n)
		.field("b", &FPCAirreg::b)
		.field("c", &FPCAirreg::c)
		.field("gap", &FPCAirreg::gap)

		.field("grid", &FPCAirreg::grid)
		.field("cuttedgrid", &FPCAirreg::cuttedgrid)

		.field("cross_sectional_mean", &FPCAirreg::cross_sectional_mean)
		.field("mean", &FPCAirreg::mean)

		.field("rawcov", &FPCAirreg::rawcov)
		.field("cov_weight", &FPCAirreg::cov_weight)
		.field("smtcov", &FPCAirreg::smtcov)
		.field("cuttedsmtcov", &FPCAirreg::cuttedsmtcov)
		.field("cuttedfitcov", &FPCAirreg::cuttedfitcov)

		.field("lambda", &FPCAirreg::lambda)
		.field("eigenfunc", &FPCAirreg::eigenfunc)

		.field("sigma", &FPCAirreg::sigma)

		;

	class_<lwls>("lwls")
		.constructor<VectorXd, VectorXd, VectorXd, VectorXd, int, int>()
		.field("x", &lwls::x)
		.field("y", &lwls::y)
		.field("w", &lwls::w)
		.field("grid", &lwls::grid)
		.field("output", &lwls::output)
		.field("new_weight", &lwls::new_weight)
		.method("optimize", &lwls::optimize)
		;

}


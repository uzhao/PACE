// must enable EIGEN2_SUPPORT otherwise we can't return MatrixXd in some case
// should be removed after they fix this problem
// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2012-May/003792.html
// public everything for debug

#define EIGEN2_SUPPORT
#define GAUSSIAN_MULT 0.3989422804014326779399460599343818684758586311649346
#define LAMBDA_THRESHOLD 0.000001

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

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

double trapz(VectorXd *x, VectorXd *y) {
	int n = x->size() - 1;
	return (x->segment(1, n) - x->segment(0, n)).cwiseProduct(y->segment(1, n) - y->segment(0, n)).sum();
}

// for mean, eigen func
// also support other Xlwls
class Lwls {
public:
	VectorXd     grid;
	VectorXd     output;
	VectorXd     new_weight;
	VectorXd     x;
	VectorXd     y;
	VectorXd     w;
	VectorXd     mean;
	VectorXd     hat;
	int          degree;
	int          drv;
	double       h;
	double       bw;
	double       minbw;
	double       maxbw;
	double       gcvv;

	Lwls(VectorXd x_, VectorXd y_, VectorXd w_, VectorXd grid_, int degree_, int drv_) : x(x_), y(y_), w(w_), grid(grid_), degree(degree_), drv(drv_) {
		output = VectorXd::Zero(grid.size());
		mean = VectorXd::Zero(grid.size());
		new_weight = VectorXd::Zero(grid.size());
		hat = VectorXd::Zero(grid.size());
		minbw = (x.tail(x.size() - 1) - x.head(x.size() - 1)).mean();
		maxbw = x[x.size() - 1] - x[0];
	};

	void optimize() {
		double best_gcv = gopt(minbw, maxbw, this, &Lwls::update);
		bw = sqrt(best_gcv * minbw);
		update(bw, true);
	};

	double update(double bw_) {
		return update(bw_, false);
	};

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
				local_x.row(i).segment(0, degree - 1) = local_x.row(i - 1).segment(1, degree - 1);
			}

			// FIXME
			// THIS HAT MATRIX IS ONLY FOR DEGREE 1
			////////////////////////////////////////////////////////////////////////////////////////////
			hat[i] = w[i] * local_x(1, 1) / (local_x(0, 0)*local_x(1, 1) - local_x(0, 1)*local_x(1, 0));
			////////////////////////////////////////////////////////////////////////////////////////////
			local_ans = local_x.jacobiSvd(ComputeThinU | ComputeThinV).solve(local_respond);
			mean[i] = local_ans[0];
			output[i] = local_ans[drv];
		}
		h = hat.mean();
		gcvv = (y - mean).array().square().sum() / (1 - 2 * h + h * h);
		if (!everything) {
			grid = true_grid;
		}
		return gcvv;
	};

	~Lwls() {
	};
};

// support other Xlwls
// smooth along column
class Mlwls {
public:
	double             bw;
	double             minbw;
	double             maxbw;
	int                degree;
	int                drv;
	double             gcvv;
	int                kernel;
	std::vector<Lwls*> lwls;
	MatrixXd           smtmat;
	MatrixXd           new_weight;

	Mlwls(VectorXd x, MatrixXd mat, MatrixXd weight, int degree_, int drv_) : 
		degree(degree_), drv(drv_) {
		smtmat = MatrixXd::Zero(mat.cols(), mat.cols());
		new_weight = MatrixXd::Zero(mat.cols(), mat.cols());
		int n = mat.cols();
		lwls = std::vector<Lwls*>(n);
		minbw = -INFINITY;
		maxbw = -INFINITY;
		for (int i = 0; i != n; i++) {
			lwls[i] = new Lwls(x, mat.col(i), weight.col(i), x, degree, drv);
			minbw = (minbw > lwls[i]->minbw) ? minbw : lwls[i]->minbw;
			maxbw = (maxbw < lwls[i]->maxbw) ? lwls[i]->maxbw : maxbw;
		}
		try {
			if (minbw > maxbw) {
				throw 0;
			}
		}
		catch (int) {
			std::cout << "The input for local weighted least square is too few!" << std::endl;
			exit(1);
		}
		optimize();
	};

	void optimize() {
		double best_gcv;
		best_gcv = gopt(minbw, maxbw, this, &Mlwls::update);
		bw = sqrt(best_gcv * minbw);
		update(bw, true);
	};

	double update(double bw_) {
		return update(bw_, false);		
	};

	double update(double bw_, bool everything) {
		bw = bw_;
		gcvv = 0;
		for (int i = 0; i != lwls.size(); i++) {
			gcvv += lwls[i]->update(bw, everything);
			if (everything) {
				smtmat.col(i) = lwls[i]->output;
				new_weight.col(i) = lwls[i]->new_weight;
			}
		}
		return gcvv;
	};

	~Mlwls() {
		for (std::vector<Lwls*>::iterator it = lwls.begin(); it != lwls.end(); it++) {
			delete *it;
		}
	};
};

// for cov mat
class Covlwls {
public:
	MatrixXd midmat;
	MatrixXd midweight;
	MatrixXd smtmat;

	Covlwls(VectorXd x, MatrixXd mat, MatrixXd weight) {
		Mlwls middle(x, mat, weight, 1, 0);
		midmat = middle.smtmat.transpose();
		midweight = middle.new_weight.transpose();

		Mlwls final(x, midmat, midweight, 1, 0);
		smtmat = (final.smtmat + final.smtmat.transpose()) / 2;
	};

	~Covlwls() {
	};
};

// diagonal elements of covariance matrix
// FIXME
//////////////////////////////////////////////////////////
class Qlwls {
public:
	VectorXd smtdiag;

	Qlwls(MatrixXd mat, MatrixXd weight, int neighbours) {
		mat = mat.rowwise().reverse();
		weight = weight.rowwise().reverse();		
	};

	~Qlwls() {
	};
};
//////////////////////////////////////////////////////////

// for partial cov mat
// can not be used until fix hat matrix in lwls and 
// figure out is bandwidth and so on
// just copy from Covlwls with new degree and drv
// class Pcovlwls {
// public:
// };

class Subject {
public:
	// speed up est raw covmat, may remove
	// VectorXd binx;
	// VectorXd biny;
	// void bin() {
	// };

	VectorXd x;
	VectorXd y;
	int      cutstart;
	int      cutend;
	double   hash;

	VectorXd pcs;

	Subject(VectorXd *x_, VectorXd *y_) : x(*x_), y(*y_){
		cutstart = Infinity;
		cutend = -Infinity;
	};

	void cutindex(double min, double max) {
		for (int i = 0; i != x.size(); i++) {
			if (x[i] >= min) {
				cutstart = i;
				break;
			}
		}
		for (int i = x.size() - 1; i != -1; i--) {
			if (x[i] <= max) {
				cutend = i;
				break;
			}
		}
		if (cutstart > cutend || cutstart == Infinity || cutend == -Infinity) {
			cutstart = Infinity;
			cutend = -Infinity;
		}
		if (cutstart != Infinity && cutend != -Infinity) {
			hash = x.segment(cutstart, cutend - cutstart + 1).array().sin().sum();
		}
		else {
			hash = Infinity;
		}
	};

	~Subject(){
	};
};

class FPCA {
public:
	List     debugvar;

	bool     regular;
	int      n;            // number of subjects
	int      b;            // number of bins
	int      c;            // number of components

	VectorXd grid;
	VectorXd cuttedgrid;
	VectorXd biny;
	VectorXd mean_weight;
	VectorXd mean;
	VectorXd cuttedmean;

	MatrixXd rawcov;
	MatrixXd cov_weight;
	MatrixXd smtcov;
	MatrixXd cuttedsmtcov;
	MatrixXd cuttedfitcov;
	
	MatrixXd *mfspline;
	Spline2d *invmfspline;
	std::vector<MatrixXd*> efspline;
	std::vector<Spline2d*> invefspline;
	std::map<double, VectorXd> invmean;
	std::map<double, MatrixXd> invef;
	std::map<double, MatrixXd> invfitcov;
	std::map<double, MatrixXd> pcsolver;

	double   cutp;
	VectorXd fve;
	double   fveth;
	int      maxk;
	VectorXd lambda;
	MatrixXd eigenfunc;

	bool     error;
	double   sigma;
	bool     new_ridge;
	double   rho;

	int      neighbours;
	int      cuttedlength;
	double   minx;
	double   maxx;
	double   cutminx;
	double   cutmaxx;
	double   gap;

	Lwls     *lwls_mean;
	Covlwls  *covlwls_cov;
	//Qlwls    *qlwls_diag;

	std::vector<Subject*> subjects;

	FPCA(List x_, List y_, double cutp_, double fveth_, int maxk_, bool error_, bool regular_) :
		cutp(cutp_), fveth(fveth_), maxk(maxk_), error(error_), regular(regular_){

		// get number of subjects
		n = x_.size();

		subjects = std::vector<Subject*>(n);

		// set number of bins
		if (n >= 400) {
			b = 50;
		}
		else if (n <= 100) {
			b = 25;
		}
		else {
			b = 2.5 * sqrt(n);
		}
		neighbours = ceil(b * cutp);

		// set subjects 
		maxx = -INFINITY;
		minx = INFINITY;

		for (int i = 0; i != n; i++) {
			VectorXd xt = x_[i];
			VectorXd yt = y_[i];
			maxx = (maxx > xt.maxCoeff()) ? maxx : xt.maxCoeff();
			minx = (minx < xt.minCoeff()) ? minx : xt.minCoeff();
			subjects[i] = new Subject(&xt, &yt);
			// subjects[i]->bin();
		}

		// to avoid overflow
		maxx *= 1 + DBL_EPSILON;

		// set grid
		gap = (maxx - minx) / b;
		grid = VectorXd::LinSpaced(b, minx + gap / 2, maxx - gap / 2);
		cuttedlength = b - 2 * neighbours;
		cuttedgrid = grid.segment(neighbours, cuttedlength);
		cutminx = cuttedgrid(0);
		cutmaxx = cuttedgrid(cuttedgrid.size() - 1);

		for (int i = 0; i != n; i++) {
			subjects[i]->cutindex(cutminx, cutmaxx);
		}

		biny = VectorXd::Zero(b);
		mean_weight = VectorXd::Zero(b);

		for (int i = 0; i != n; i++) {
			for (int j = 0; j != subjects[i]->x.size(); j++) {
				int index = (subjects[i]->x[j] - minx) / gap;
				(mean_weight[index])++;
				biny[index] += (subjects[i]->y[j] - biny[index]) / mean_weight[index];
			}
		}

		// get mean
		lwls_mean = new Lwls(grid, biny, mean_weight, grid, 1, 0);
		lwls_mean->optimize();
		mean = lwls_mean->output;
		cuttedmean = mean.segment(neighbours, cuttedlength);

		// get rawcov
		// no binning for estimating raw cov matrix
		cov_weight = MatrixXd::Zero(b, b);
		rawcov = MatrixXd::Zero(b, b);
		for (int i = 0; i != n; i++) {
			int m = subjects[i]->x.size();
			for (int j = 0; j != m; j++) {
				int index_j = (subjects[i]->x[j] - minx) / gap;
				for (int k = j; k != m; k++) {
					int index_k = (subjects[i]->x[k] - minx) / gap;
					cov_weight(index_j, index_k)++;
					rawcov(index_j, index_k) += ((subjects[i]->y[j] - mean[index_j]) * (subjects[i]->y[k] - mean[index_k]) - rawcov(index_j, index_k)) / cov_weight(index_j, index_k);
				}
			}
		}

		cov_weight = cov_weight + cov_weight.transpose();
		// cov_weight.diagonal() = cov_weight.diagonal() / 2;
		VectorXd oricov_weight = cov_weight.diagonal().segment(neighbours, cuttedlength) / 2;
		cov_weight.diagonal() *= 0;
		rawcov = rawcov + rawcov.transpose();
		rawcov.diagonal() = rawcov.diagonal() / 2;

		// get smtcov
		covlwls_cov = new Covlwls(grid, rawcov, cov_weight);
		smtcov = covlwls_cov->smtmat;

		/*
		  TODO: Weight? Before eigen or after? Number of nbs?
		  fix smtcov diag with qdiag
		  */
		// qlwls_diag = new Qlwls(rawcov, cov_weight, neighbours);
		// smtcov.diagonal() = qlwls_diag->smtdiag;

		cuttedsmtcov = smtcov.block(neighbours, neighbours, cuttedlength, cuttedlength);

		// eigen decomposition target cov matrix
		SelfAdjointEigenSolver<MatrixXd> eigendcp(cuttedsmtcov);
		lambda = eigendcp.eigenvalues().reverse();
		int posk = 0;
		double lambdasum = (lambda.sum() + lambda.cwiseAbs().sum()) / 2;
		for (int i = 0; i != lambda.size(); i++) {
			if (lambda[i] / lambdasum > LAMBDA_THRESHOLD) {
				posk++;
			}
		}
		lambda = lambda.head(posk) * gap;
		eigenfunc = eigendcp.eigenvectors().rowwise().reverse().leftCols(posk) / sqrt(gap);

		// choose number of components (at least 2)
		c = 2;
		fve = lambda / lambda.sum();
		for (int i = 1; i != fve.size(); i++) {
			fve[i] += fve[i - 1];
			if (fve[i] < fveth) {
				c = i + 2;
			}
		}
		c = c > maxk ? maxk : c;

		// fit cov matrix for group
		cuttedfitcov = MatrixXd::Zero(cuttedlength, cuttedlength);
		for (int i = 0; i != c; i++) {
			cuttedfitcov += lambda[i] * (eigenfunc.col(i) * eigenfunc.col(i).transpose());
		}

		// est sigma
		Lwls covdiag = Lwls(cuttedgrid, rawcov.diagonal().segment(neighbours, cuttedlength), oricov_weight, cuttedgrid, 1, 0);
		covdiag.optimize();
		VectorXd vare = covdiag.output - cuttedsmtcov.diagonal();
		sigma = trapz(&cuttedgrid, &vare) / (cuttedgrid[cuttedlength - 1] - cuttedgrid[0]);
		sigma = sigma > 0 ? sigma : 0;

		//// construct spline for mean function
		mfspline = new MatrixXd;
		mfspline->resize(2, cuttedlength);
		mfspline->row(0) = cuttedgrid;
		mfspline->row(1) = cuttedmean;
		invmfspline = new Spline2d;
		(*invmfspline) = SplineFitting<Spline2d>::Interpolate((*mfspline), 3);

		// construct spline for eigen function
		efspline.resize(c);
		invefspline.resize(c);
		for (int i = 0; i != c; i++) {
			efspline[i] = new MatrixXd;
			efspline[i]->resize(2, cuttedlength);
			efspline[i]->row(0) = cuttedgrid;
			efspline[i]->row(1) = eigenfunc.col(i).transpose();
			invefspline[i] = new Spline2d;
			*(invefspline[i]) = SplineFitting<Spline2d>::Interpolate(*(efspline[i]), 3);
		}

		// fit mean, 
		// fit eigen function, 
		// fit cov matrix,
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			if ((*it)->hash == Infinity) {
				continue;
			}

			if (invmean.find((*it)->hash) == invmean.end()) {
				invmean[(*it)->hash] = fitmean(*it);
			//}
			//if (invef.find((*it)->hash) == invef.end()) {
				invef[(*it)->hash] = fitef(*it);
			//}
			//if (invfitcov.find((*it)->hash) == invfitcov.end()) {
				invfitcov[(*it)->hash] = fitcov(&invef[(*it)->hash]);
			//}
			//if (pcsolver.find((*it)->hash) == pcsolver.end()) {
				MatrixXd capsigma = invfitcov[(*it)->hash];
				capsigma.diagonal() = capsigma.diagonal().array() + sigma;
				pcsolver[(*it)->hash] = lambda.head(c).asDiagonal() * invef[(*it)->hash] * capsigma.inverse();
			}
		}

		// update ridge
		// not implement yet
		// if (new_ridge) {
		//	 update_ridge();
		// }

		// estimate pc scroe
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			if ((*it)->hash == Infinity) {
				continue;
			}

			(*it)->pcs = (pcsolver[(*it)->hash] *
				          ((*it)->y.segment((*it)->cutstart, (*it)->cutend - (*it)->cutstart + 1) - 
						   invmean[(*it)->hash])).array();
			debugvar.push_back((*it)->pcs);
		}


	};

	~FPCA() {
		for (std::vector<Subject*>::iterator it = subjects.begin(); it != subjects.end(); it++) {
			delete *it;
		}
		//for (std::vector<MatrixXd*>::iterator it = efspline.begin(); it != efspline.end(); it++) {
		//	delete *it;
		//}
		//for (std::vector<Spline2d*>::iterator it = invefspline.begin(); it != invefspline.end(); it++) {
		//	delete *it;
		//}
		//for (std::map<double, MatrixXd*>::iterator it = invfitcov.begin(); it != invfitcov.end(); it++) {
		//	delete it->second;
		//}
		delete mfspline;
		delete invmfspline;
		delete lwls_mean;
		delete covlwls_cov;
		//delete qlwls_diag;
	};

	void update_ridge() {
	};

	VectorXd fitline(Subject *it, Spline2d *s) {
		int size = it->cutend - it->cutstart + 1;
		VectorXd ans;
		ans.setZero(size);
		for (int i = 0; i != size; i++) {
			double w = it->x(it->cutstart + i);
			ans(i) = (*s)((w - cutminx) / (cutmaxx - cutminx))(1);
		}
		return ans;
	};

	VectorXd fitmean(Subject *it) {
		return fitline(it, invmfspline);
	};

	MatrixXd fitef(Subject *it) {
		MatrixXd ans = MatrixXd::Zero(c, it->cutend - it->cutstart + 1);
		for (int i = 0; i != c; i++) {
			ans.row(i) = fitline(it, invefspline[i]);
		}
		return ans;
	};

	MatrixXd fitcov(MatrixXd *it) {
		MatrixXd ans = MatrixXd::Zero(it->cols(), it->cols());
		for (int i = 0; i != it->rows(); i++) {
			ans += it->row(i).transpose() * it->row(i);
		}
		return ans;
	};

	Subject *index(int i) {
		return subjects[i - 1];
	};


};

RCPP_MODULE(PACE){
	using namespace Rcpp;

	class_<FPCA>("FPCA")
		.constructor<List, List, double, double, int, bool, bool>()

		.field("n", &FPCA::n)
		.field("b", &FPCA::b)
		.field("c", &FPCA::c)

		.field("grid", &FPCA::grid)
		.field("mean", &FPCA::mean)

		.field("rawcov", &FPCA::rawcov)
		.field("cov_weight", &FPCA::cov_weight)
		.field("smtcov", &FPCA::smtcov)
		.field("cuttedsmtcov", &FPCA::cuttedsmtcov)

		.field("cuttedfitcov", &FPCA::cuttedfitcov)

		.field("lambda", &FPCA::lambda)
		.field("eigenfunc", &FPCA::eigenfunc)
		.field("fve", &FPCA::fve)

		.field("sigma", &FPCA::sigma)

		.field("biny", &FPCA::biny)
		.field("gap", &FPCA::gap)

		.field("debugvar", &FPCA::debugvar)

		.method("index", &FPCA::index)
		;
}


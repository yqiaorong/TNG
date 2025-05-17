import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.optimize import fmin, minimize
from math import isnan

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'font.family' : 'serif', 'font.size' : 12,
                 'text.latex.preamble' : r"\usepackage{amsmath}",
                 'mathtext.fontset' : 'dejavuserif', 'xtick.major.pad' : 2,
                 'ytick.major.pad' : 2, 'xtick.major.size' : 6,
                 'ytick.major.size' : 6, 'xtick.minor.size' : 3,
                 'ytick.minor.size' : 3, 'axes.linewidth' : 2, 'axes.labelpad' : 1})
plt.rcParams["text.usetex"] = False

class splashback_properties:
    def __init__(self, tag, x, y, mean, Ob=0.0486, Om=0.3089, R200mean=None, INNER='EINASTO', \
                 OUTER='MATTER', Comparison=False, Norm=True, smooth=False):


        self.tag = tag
        self.Comparison = Comparison
        print('  -Computing {0} splashback properties'.format(self.tag))
        # Store non-zero profiles all work will be done with -- assumes these are NOT logged!!!
        idx = np.where(y != 0.0)[0]
        self.x = x[idx]
        self.y = y[idx]

        # Store MEAN density of Universe, compute R200mean
        self.rho_mean = mean
        self.Ob       = Ob
        self.Om       = Om
        if R200mean is None:
            self.R200mean = self.compute_R200_mean()
        else:
            self.R200mean = R200mean

        # Set inner and outer profiles to be used
        self.inner_pro = INNER
        self.outer_pro = OUTER


        # Compute unsmoothed and smoothed numerical derivatives
        self.profile_num_deriv = self.compute_numerical_deriv(self.x, self.y)
        self.global_min_x = self.x[np.argmin(self.profile_num_deriv)]

        if smooth:
            self.do_smoothing()
            self.smooth_num_deriv = self.compute_numerical_deriv(self.x, self.y, SG=True)

        return


    def do_smoothing(self):
        """
        Smooth the data with Savitsky-Golay, find Rspl
        """

        wind = 15
        poly = 4

        radius  = self.x
        density = self.y


        # savgol needs evenly spaced data points - for smoothing in log space, interpolate between points
        # Check if points are even
        logspaced = ( (np.log10(radius[2])-np.log10(radius[1])) == (np.log10(radius[1])-np.log10(radius[0])) )
        if not logspaced:
            interp_to_logspace = interp1d(np.log10(radius),np.log10(density))
            radius = np.logspace(np.log10(radius[0]),np.log10(radius[-1]),int(len(radius)/2))
            density = 10.**interp_to_logspace(np.log10(radius))

        # Spacing between points
        dlt = np.log10(radius[1])-np.log10(radius[0])

        # Logspaced smoothing and derivatives
        log_radius = np.log10(radius)
        log_density = np.log10(density)
        smoothed_density = 10**savgol_filter(log_density,wind,poly,delta=dlt)
        smoothed_slope = savgol_filter(log_density,wind,poly,delta=dlt,deriv=1)
        deriv_slope = savgol_filter(log_density,wind,poly,delta=dlt,deriv=2)
        deriv_slope[1:] = savgol_filter(log_density[1:],wind,poly,delta=2,deriv=2)

        # Steepest slope is where second derivative = 0.  Interpolate between points to get to 0
        smoothedmin = np.argmin(smoothed_slope[3:-2])+3
        interpderiv = interp1d(log_radius,deriv_slope)
        interp_radius = np.linspace(log_radius[smoothedmin-1],log_radius[smoothedmin+1])
        interp_deriv_slope = interpderiv(interp_radius)
        min_slope_ind = np.argmin(np.abs(interp_deriv_slope))

        log_Rspl = interp_radius[min_slope_ind]
        Rspl = 10**log_Rspl

        self.smooth_radius  = radius
        self.smooth_density = smoothed_density
        self.smooth_deriv   = smoothed_slope
        self.smooth_Rsp     = Rspl
        self.smoothed       = True

        return


    def compute_R200_mean(self):
        """
        Use the density profile to compute R200 mean. NOTE: this is not a
        fantastic way to compute this and only technically works for with
        the total matter profile, hence the option to supply it in init
        """

        Rint = interp1d(self.y / self.rho_mean, self.x)
        return Rint(200.0)

    def compute_numerical_deriv(self, x, y, log=True, SG=False, order=4):
        """
        This function computes the numerical derivative assuming first order
        finite different scheme -- NOTE: Computation in log10 space
        """

        if SG: y = savgol_filter(y, 5, 1)

        if log:
            lgR = np.log10(x)
            lgP = np.log10(y)
        else:
            lgR = x
            lgP = y

        dlgP_dlgR       = np.zeros(len(lgR), dtype=np.float32)

        if order==1:
            dlgP_dlgR[1:-1] = (lgP[2:] - lgP[:-2]) / (lgR[2:] - lgR[:-2])
            dlgP_dlgR[0]    = (lgP[1] - lgP[0]) / (lgR[1] - lgR[0])
            dlgP_dlgR[-1]   = (lgP[-1] - lgP[-2]) / (lgR[-1] - lgR[-2])

        elif order==4:
            dlgP_dlgR[2:-2] = ( 1./12.*lgP[0:-4] - 2./3.*lgP[1:-3] + 2./3.*lgP[3:-1] - 1./12.*lgP[4:] ) / (lgR[3:-1] - lgR[2:-2])
            dlgP_dlgR[0]    = (lgP[1] - lgP[0]) / (lgR[1] - lgR[0])
            dlgP_dlgR[1]    = (lgP[2] - lgP[0]) / (lgR[2] - lgR[0])
            dlgP_dlgR[-1]   = (lgP[-1] - lgP[-2]) / (lgR[-1] - lgR[-2])
            dlgP_dlgR[-2]   = (lgP[-3] - lgP[-1]) / (lgR[-3] - lgR[-1])
            #dlgP_dlgR[0:2] += dlgP_dlgR[2]
            #dlgP_dlgR[-2:] += dlgP_dlgR[-3]

        del lgP, lgR

        return dlgP_dlgR

    def fit_DK14(self, p0=None, CHECK=True, selection=None, peak_height=None, accretion_rate=None):

        self.selection = selection
        self.peak_height = peak_height
        self.accretion_rate = accretion_rate


        """
                p0 = [ 2.01621038e+04,  1.00878323e-02,  5.24792318e+00,  6.60942272e-01, 6.59144753e-01,  4.36929551e+00,  5.49451986e+04,  1.75506745e-01, -9.99990122e+00,  1.07008381e-08,  4.59741827e+00]
        Fit the DK14 profile to the provided profile
        """
        CHECK = self.Comparison
        #--- Initial guess and limits
        if self.inner_pro in ['EINASTO', 'BETA']:
            free_params = [0,1,3,-2,-1]
            if p0 is None:
                # rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e
                p0 = [1.0e4 * self.rho_mean, 0.5 * self.R200mean, 1.0, 1.0 * self.R200mean, 2.0, 4.0, 0.0, 0.0]
                # p0 for fitting profile r^2
                #p0 = [288601.87942122505, 0.15125536617648788, 1., 1.102902248833148, 2., 4., 0.4928676055716866, 1.8946642725751321]
                #p0 = [286075.60010078945 * (self.rho_mean / 39.339), 0.13272954836296388 * self.R200mean, 0.2262245729158164, 1.5999129096532894 * self.R200mean, 7.959999489473648, 9.999999999999998, 0.45103488531622626, 2.138749804789657]
                if self.selection is not None:
                    p0 = [p0[i] for i in free_params]
                    #p0 = [6.10520502e+05, 9.34768805e-02, 1.45816596e+00, 4.25129755e-01, 1.59129287e+00]
                    #p0 = [self.y[0]/np.e**2, .01 * self.R200mean, self.R200mean, 0.5, 1.5]
                    print(p0,'\n\n')
#            limits = ([0.0, 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0],
##                      [1.0e2*p0[0], 1.0e2 * self.R200mean, 1.0e2*self.R200mean, 10., 10.])
#                      [np.inf, np.inf, 10.0, np.inf, 10.0, 10.0, 10.0, 10.0])
##                      [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])
            limits = (#[0.0, 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0],
                      #[np.inf, np.inf, 10.0, np.inf, 10.0, 20.0, 10.0, 10.0])
                      [np.min(self.y), np.min(self.x), 0.,  np.min(self.x), 0., 0.,  0.,  0.],
                      [10.*np.max(self.y), np.max(self.x), 10., np.max(self.x), 10.,20., 10., 10.])


            if self.selection is not None:
                    limits = [[limits[0][i] for i in free_params], [limits[1][i] for i in free_params]]
        elif self.inner_pro in ['VIK06']:
            if p0 is None:
                # rho_s, r_s, alpha, b, r_x, e, r_t, beta, gamma, b_e, s_e
                p0 = [1000.0 * self.rho_mean, 0.1 * self.R200mean, 2.0, 2.0, \
                      0.8 * self.R200mean, 2.5, 1.0 * self.R200mean, 2.0, 4.0, 0.0, 0.0]
            limits = ([0.1, 0.01 * self.R200mean, 0.0, 0.0, 0.01 * self.R200mean, 0.0,
                       0.0, -10.0, -10.0, -10.0, -10.0],
                      [1.0e6, np.inf, 10.0, np.inf, np.inf, 10.0,
                       np.inf, 10.0, 5.0, 10.0, 10.0])

        # Actual fit
        p0 = np.clip(p0, limits[0], limits[1])
        fit = least_squares(self.DK14_function, p0, args=(self.inner_pro, self.outer_pro, 'FIT'), \
                            bounds=limits, method='trf')

        if self.selection is not None:
            fit_count = 1
            while fit.success == False:
                print('Try',fit_count)
                p0 = fit.x
                fit = least_squares(self.DK14_function, p0, args=(self.inner_pro, self.outer_pro, 'FIT'), \
                                    bounds=limits, method='trf')
                if fit_count == 100:
                    print('100 fits!\n')
                    break
                fit_count += 1

        #print(fit)
        #print(p0)
        # Store parameters
        if self.inner_pro in ['EINASTO', 'BETA']:
            if self.selection is None:
                self.rho_s, self.r_s, self.alpha, self.r_t, self.beta, self.gamma, \
                    self.b_e, self.s_e = fit.x
            else:
                self.rho_s, self.r_s, self.r_t, self.b_e, self.s_e = fit.x
        elif self.inner_pro in ['VIK06']:
            self.rho_s, self.r_s, self.alpha, self.b, self.r_x, self.e, \
                self.r_t, self.beta, self.gamma, self.b_e, self.s_e = fit.x 

        self.fit_success = fit.success
        self.fit_params  = fit.x


        self.fit = fit
        print("     message:",fit.message)
        print("        nfev:",fit.nfev)
        print("        njev:",fit.njev)
        print("  optimality:",fit.optimality)
        print("      status:",fit.status)
        print("     success:",fit.success)
        print("           x:",fit.x)
        print()


        # Check fit -- if required
        if CHECK:
            model = self.DK14_function(fit.x, self.inner_pro, self.outer_pro, 'MODEL')

            fig = plt.figure(figsize=(5,5))
            ax  = fig.add_subplot(111)
            pl1 = ax.loglog(self.x, self.y / self.rho_mean, 's', c='darkgrey', ms=4, \
                            label='Data', zorder=-4)
            pl2 = ax.loglog(self.x, model / self.rho_mean, '-', c='crimson', lw=2, \
                            label='Fit', zorder=-3)
            if hasattr(self, 'smooth_density'):
                pl3 = ax.loglog(self.smooth_radius, self.smooth_density / self.rho_mean, '-', c='black', lw=2, \
                                label='Smooth')
            if np.min(self.x) < 1.:
                ax.set_xlim(.009,5.5)
                ax.set_xlabel(r'$r/r_{200m}$', fontsize=18, family='serif')
            else:
                ax.set_xlim(1.0e1, 1.0e4)
                ax.set_xlabel(r'$r$  $[\mathrm{kpc}]$', fontsize=18, family='serif')
            ax.set_ylim(0.1, 5.0e6)
            ax.set_ylabel(r'$\rho\,/\,\rho_{\mathrm{mean}}$', fontsize=18, family='serif')
            ax.legend(loc=1, fontsize=10)
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            fig.tight_layout()
            fig.savefig('images/{0}_density_comparison.png'.format(self.tag), dpi=128)
            plt.close(fig)
        return


    def DK14_function(self, p0, INNER, OUTER, MODE, Ob=0.0486, Om=0.3089):
        """
        Compute the DK14 model based on parameters supplied
        """

        # Inner
        if INNER == 'EINASTO':
            if len(p0) == 5:
                rho_s, r_s, r_t, b_e, s_e = p0
                if self.peak_height is None:
                    print('Need peak_height for fitting with selection option!')
                    quit()
                alpha = 0.155 + 0.0095 * self.peak_height**2
                if  self.selection == 'mass':
                    #r_t = (1.9 - 0.18*self.peak_height) * self.R200mean
                    beta = 4
                    gamma = 8
                elif self.selection == 'accretion':
                    #if self.accretion_rate is None:
                    #    print('Need accretion_rate for selection option "accretion"')
                    #    quit()
                    #r_t = (0.62 + 1.18 * np.exp(-self.accretion_rate)) * self.R200mean
                    beta = 6
                    gamma = 4
                self.alpha = alpha
                #self.r_t = r_t
                self.beta = beta
                self.gamma = gamma
            else:
                rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e = p0
            rho_inner = rho_s * np.exp(-(2.0 / alpha) * ((self.x / r_s) ** alpha - 1.0))
        elif INNER == 'BETA':
            rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e = p0
            rho_inner = rho_s / (1.0 + (self.x / r_s) ** 2.0) ** (1.5 * alpha)
        elif INNER == 'VIK06':
            rho_s, r_s, alpha, b, r_x, e, r_t, beta, gamma, b_e, s_e = p0
            rho_inner = rho_s * ((self.x / r_s) ** (-alpha * 0.5) \
                                / (1.0 + (self.x / r_s) ** 2.0) ** (3.0 * b / 2.0 - alpha / 4.0)) \
                * (1.0 / ((1.0 + (self.x / r_x) ** 3.0) ** (e / 6.0)))
        else:
            print('ERROR:\n---> {0} inner profile not implemented yet!\nEXITING'.format(INNER))
            quit()

        # Transition
        rho_trans = (1.0 + (self.x / r_t) ** beta) ** (-gamma / beta)

        # Outer
        if OUTER == 'MATTER':
            rho_outer = self.rho_mean * (b_e * (self.x / (5.0 * self.R200mean)) ** -s_e + 1.0)
        elif OUTER == 'BARYONIC':
            rho_outer = (self.Ob / self.Om) * self.rho_mean * \
                        (b_e * (self.x / (5.0 * self.R200mean)) ** -s_e + 1.0)
        else:
            print('ERROR:\n---> {0} outer profile not implemented yet!\nEXITING'.format(INNER))
            quit()

        # Compute total model
        model = rho_inner * rho_trans + rho_outer

        # Return required output -- fit returns r^2 * model to reduce dynamic range
        if MODE == 'FIT':
            ydx = np.where(self.y <= 0)[0]
            mdx = np.where(model <= 0)[0]
            if len(ydx) >= 1 or len(mdx) >= 1: return np.zeros(len(model)) + np.inf
            return np.log10(model / self.y)
            #return ( (model - self.y) * self.x**2 )
        elif MODE == 'MODEL':
            return model
        elif MODE == 'NORM':
            return model*self.density_normalization
        return


    def fit_DK14slope(self, p0=None, selection=None, peak_height=None, accretion_rate=None):
        '''
        Get fit parameters by fitting to the derivative of the density profile
        '''

        self.selection = selection
        self.peak_height = peak_height
        self.accretion_rate = accretion_rate


        #--- Initial guess and limits
        if self.inner_pro in ['EINASTO', 'BETA']:
            if p0 is None:
                # rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e
                #p0     = [1.0e4 * self.rho_mean, 0.5 * self.R200mean, 1.0, 1.0 * self.R200mean, 2.0, 4.0, 0.0, 0.0]
                p0 = [286075.60010078945 * (self.rho_mean / 39.339), 0.13272954836296388 * self.R200mean, 0.2262245729158164, 1.5999129096532894 * self.R200mean, 7.959999489473648, 9.999999999999998, 0.45103488531622626, 2.138749804789657]
            limits = (#[0.0, 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0],
                      #[np.inf, np.inf, 10.0, np.inf, 10.0, 20.0, 10.0, 10.0])
                      [np.min(self.y), np.min(self.x), 0.,  np.min(self.x), 0., 0.,  0.,  0.],
                      [10.*np.max(self.y), np.max(self.x), 10., np.max(self.x), 10.,20., 10., 10.])
            if self.selection is not None:
                # rho_s, r_s, r_t, b_e, s_e
                p0_inds = [0,1,3,6,7]
                p0 = [p0[i] for i in p0_inds]
                limits = [ [limits[0][i] for i in p0_inds], [limits[1][i] for i in p0_inds] ]
        elif self.inner_pro in ['VIK06']:
            # USE THIS FOR MASS STACKING
            if p0 is None:
                # rho_s, r_s, alpha, b, r_x, e, r_t, beta, gamma, b_e, s_e
                p0     = [1.0e4 * self.rho_mean, 0.1 * self.R200mean, 4.0, 1.2, \
                          0.8 * self.R200mean, 2.5, 1.0 * self.R200mean, 2.0, 4.0, 0.0, 0.0]

            limits = ([0.1, 0.01 * self.R200mean, 0.0, 0.0, 0.01 * self.R200mean, 0.0,
                       0.0, -10.0, -10.0, -10.0, -10.0],
                      [1.0e6*self.rho_mean, 1.0e3*self.R200mean, 20.0, 20.0, 1.0e3*self.R200mean, 20.0,
                       1.0e3*self.R200mean, 20.0, 20.0, 20.0, 20.0])

            # USE THIS FOR ACCRETION STACKING
            #if p0 is None:
            #    # rho_s, r_s, alpha, b, r_x, e, r_t, beta, gamma, b_e, s_e
            #    #p0 = [1.0e4 * self.rho_mean, 0.1 * self.R200mean, 4.0, 1.2, \
            #    #          0.8 * self.R200mean, 2.5, 1.0 * self.R200mean, 2.0, 4.0, 0.0, 0.0]
            #    p0 = [1.0e3 * self.rho_mean, 0.1 * self.R200mean, 4.0, 1.2, \
            #              0.8 * self.R200mean, 2.5, 1.0 * self.R200mean, 2.0, 4.0, 0.0, 0.0]
            #limits = (#[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],#[0.1, 0.01 * self.R200mean, 0.0, 0.0, 0.01 * self.R200mean, 0.0,
            #           #0.0, -0.0, -0.0, -0.0, -0.0],
            #          [np.min(self.y), np.min(self.x), 0., 0., np.min(self.x), 0., np.min(self.x), 0., 0., 0., 0.],
            #          #[1.0e6*self.rho_mean, 1.0e3*self.R200mean, 20.0, 20.0, 1.0e3*self.R200mean, 20.0,
            #          # 1.0e3*self.R200mean, 20.0, 100.0, 20.0, 20.0])
            #          [10.*np.max(self.y), np.max(self.x), 10., 10., np.max(self.x), 10., np.max(self.x), 10., 20., 10., 10.])
            #          #[np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])
            #          # Toggle Beta between <10 and <20 if second thin dip causes problems
            #          # Use gamma <10 for 1250 resolution and <20 for 2500

        # Actual fit
        try:
            fit = least_squares(self.compute_analytic_deriv, p0, args=('FIT','dummy'), bounds=limits, method='trf')
        except:
            print('\n\nFAILING LEAST SQUARES',flush=True)
            print(p0)
            print(self.x)
            print(self.y)
            try:
                p0_o = [286075.60010078945 * (self.rho_mean / 39.339), 0.13272954836296388 * self.R200mean, 0.2262245729158164, 1.5999129096532894 * self.R200mean, 7.959999489473648, 9.999999999999998, 0.45103488531622626, 2.138749804789657]
                fit = least_squares(self.compute_analytic_deriv, p0_o, args=('FIT','dummy'), bounds=limits, method='trf')
                print('SUCCESSFUL RETRY')
            except:
                fit = least_squares(self.compute_analytic_deriv, p0, args=('FIT','dummy'), bounds=limits, method='trf')

        # Store parameters
        if self.inner_pro in ['EINASTO', 'BETA']:
            if self.selection is None:
                self.rho_s, self.r_s, self.alpha, self.r_t, self.beta, self.gamma, \
                    self.b_e, self.s_e = fit.x
            else:
                self.rho_s, self.r_s, self.r_t, self.b_e, self.s_e = fit.x

        elif self.inner_pro in ['VIK06']:
            self.rho_s, self.r_s, self.alpha, self.b, self.r_x, self.e, \
                self.r_t, self.beta, self.gamma, self.b_e, self.s_e = fit.x

        self.fit_success = fit.success
        self.fit_params  = fit.x
        self.fit = fit
        print("     message:",fit.message)
        print("        nfev:",fit.nfev)
        print("        njev:",fit.njev)
        print("  optimality:",fit.optimality)
        print("      status:",fit.status)
        print("     success:",fit.success)
        print("           x:",fit.x)
        print()


        #if self.inner_pro == 'VIK06':
            #print([self.fit_params[0]/self.rho_mean, self.fit_params[1]/self.R200mean, self.fit_params[2], self.fit_params[3], \
            #       self.fit_params[4]/self.R200mean, self.fit_params[5], self.fit_params[6], self.fit_params[7], self.fit_params[8], \
            #       self.fit_params[9], self.fit_params[10]])
            #print(fit)
        #print('\n')


        # normalize the density profile
        def normalize_density(x):
            return np.log10(x*self.DK14_function(p0=self.fit_params,INNER=self.inner_pro, OUTER=self.outer_pro, MODE='MODEL') / self.y)

        mid_index = int(len(self.y)/2)
        while isnan(self.y[mid_index]):
            mid_index += 1
        try:
            #print('\n\n')
            #print('\n\n', flush=True)
            #print(mid_index, flush=True)
            #print(self.y[mid_index], flush=True)
            #print(self.x, flush=True)
            #print(self.y, flush=True)
            #print(self.fit_params, flush=True)
            #print(self.inner_pro, flush=True)
            #print(self.outer_pro, flush=True)
            normalize_fit = least_squares(normalize_density, \
        	                    self.y[mid_index]/self.DK14_function(p0=self.fit_params,INNER=self.inner_pro, OUTER=self.outer_pro, MODE='MODEL')[mid_index]) 
        except:
            print('\n\n', flush=True)
            print('FAILING LEAST SQUARES')
            print(mid_index)
            print(self.y[mid_index], flush=True)
            print(self.x, flush=True)
            print(self.y, flush=True)
            print(self.fit_params, flush=True)
            print(self.inner_pro, flush=True)
            print(self.outer_pro, flush=True)
            normalize_fit = least_squares(normalize_density, \
                                    self.y[mid_index]/self.DK14_function(p0=self.fit_params,INNER=self.inner_pro, OUTER=self.outer_pro, MODE='MODEL')[mid_index]) 

        self.density_normalization = normalize_fit.x[0]



        if self.Comparison:
            model = self.DK14_function(fit.x, self.inner_pro, self.outer_pro, 'NORM')

            fig = plt.figure(figsize=(5,5))
            ax  = fig.add_subplot(111)
            pl1 = ax.loglog(self.x, self.y / self.rho_mean, 's', c='darkgrey', ms=4, \
                            label='Data', zorder=-4)
            pl2 = ax.loglog(self.x, model / self.rho_mean, '-', c='crimson', lw=2, \
                            label='Fit', zorder=-3)
            if hasattr(self, 'smooth_density'):
                pl3 = ax.loglog(self.smooth_radius, self.smooth_density / self.rho_mean, '-', c='black', lw=2, \
                                label='Smooth')
            if np.min(self.x) < 1.:
                ax.set_xlim(.009,5.5)
                ax.set_xlabel(r'$r/r_{200m}$', fontsize=18, family='serif')
            else:
                ax.set_xlim(1.0e1, 1.0e4)
                ax.set_xlabel(r'$r$  $[\mathrm{kpc}]$', fontsize=18, family='serif')
            ax.set_ylim(0.1, 5.0e6)
            ax.set_ylabel(r'$\rho\,/\,\rho_{\mathrm{mean}}$', fontsize=18, family='serif')
            ax.legend(loc=1, fontsize=10)
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            fig.tight_layout()
            fig.savefig('images/{0}_density_comparison.png'.format(self.tag), dpi=128)
            plt.close(fig)



        return




    def DK14_derivatives(self):
        """
        Compute the numerical derivative of the fit function and the analytic derivative
        """

        # Numerical derivative of the fit
        if self.inner_pro in ['EINASTO', 'BETA']:
            p0 = [self.rho_s, self.r_s, self.alpha, self.r_t, self.beta, self.gamma, \
                  self.b_e, self.s_e]
        elif self.inner_pro in ['VIK06']:
            p0 = [self.rho_s, self.r_s, self.alpha, self.b, self.r_x, self.e, \
                  self.r_t, self.beta, self.gamma, self.b_e, self.s_e]
        model = self.DK14_function(p0, self.inner_pro, self.outer_pro, 'MODEL')

        self.fit_num_deriv = self.compute_numerical_deriv(self.x, model)

        # Analytic derivative via fit parameters
        self.analytic_deriv = self.compute_analytic_deriv()
        return

    def compute_analytic_deriv(self, p0=None, MODE='MODEL',dummy='dummy'):
        """
        Use the fit parameters to compute the analytic derivative
        """

        if p0 is None:
            if self.inner_pro in ['EINASTO', 'BETA']:
                # rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e
                p0 = [self.rho_s, self.r_s, self.alpha, self.r_t, self.beta, self.gamma, self.b_e, self.s_e]
            elif self.inner_pro == 'VIK06':
                p0 = [self.rho_s, self.r_s, self.alpha, self.b, self.r_x, self.e, self.r_t, self.beta, self.gamma, self.b_e, self.s_e]
            else:
                print('ERROR:\n---> {0} inner profile not implemented yet!\nEXITING'.format(INNER))
                quit()

        # Original functions
        if self.inner_pro == 'EINASTO':
            if len(p0) == 5:
                rho_s, r_s, r_t, b_e, s_e = p0
                if self.peak_height is None:
                    print('Need peak_height for fitting with selection option!')
                    quit()
                alpha = 0.155 + 0.0095 * self.peak_height**2
                if  self.selection == 'mass':
                    #r_t = (1.9 - 0.18*self.peak_height) * self.R200mean
                    beta = 4
                    gamma = 8
                elif self.selection == 'accretion':
                    #if self.accretion_rate is None:
                    #    print('Need accretion_rate for selection option "accretion"')
                    #    quit()
                    #r_t = (0.62 + 1.18 * np.exp(-self.accretion_rate)) * self.R200mean
                    beta = 6
                    gamma = 4
                self.alpha = alpha
                #self.r_t = r_t
                self.beta = beta
                self.gamma = gamma
            else:
                rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e = p0

            rho_inner = rho_s * np.exp(-2.0 / alpha \
                                            * ((self.x / r_s) ** alpha - 1.0))
        elif self.inner_pro == 'BETA':
            rho_s, r_s, alpha, r_t, beta, gamma, b_e, s_e = p0
            rho_inner = rho_s / (1.0 + (self.x / r_s) ** 2.0) ** (1.5 * alpha)
        elif self.inner_pro == 'VIK06':
            rho_s, r_s, alpha, b, r_x, e, r_t, beta, gamma, b_e, s_e = p0
            rho_inner = rho_s * ((self.x / r_s) ** (-alpha * 0.5) \
                                / (1.0 + (self.x / r_s) ** 2.0) \
                                      ** (3.0 * b / 2.0 - alpha / 4.0)) \
                * (1.0 / ((1.0 + (self.x / r_x) ** 3.0) ** (e / 6.0)))

        ftrans    = (1.0 + (self.x / r_t) ** beta) ** (-gamma / beta)

        if self.outer_pro == 'MATTER':
            rho_outer = self.rho_mean * \
                        (b_e * (self.x / (5.0 * self.R200mean)) ** -s_e + 1.0)
        elif self.outer_pro == 'BARYONIC':
            rho_outer = (self.Ob / self.Om) * self.rho_mean * \
                        (b_e * (self.x / (5.0 * self.R200mean)) ** -s_e + 1.0)

        profile = rho_inner * ftrans + rho_outer

        # Gradients
        if self.inner_pro == 'EINASTO':
            d_rho_inner = -2.0 / r_s * (self.x / r_s) ** (alpha - 1.0) * rho_inner
        elif self.inner_pro == 'BETA':
            d_rho_inner = -3.0 * rho_s * alpha * self.x \
                          * ((1.0 + (self.x / r_s) ** 2.0) ** (-1.5 * alpha - 1.0)) \
                          / r_s ** 2.0
        elif self.inner_pro == 'VIK06':
            term1   = rho_s * (self.x/r_s) ** (-alpha * 0.5)
            term2   = ( 1.0 + (self.x/r_s)**2.0 ) ** (-3.0*b/2.0 + alpha/4.0)
            term3   = ( 1.0 + (self.x/r_x)**3.0 ) ** (-e/6.0)
            d_term1 = -alpha/2.0 * (self.x/r_s) ** (-alpha*0.5 - 1.0) * rho_s / r_s
            d_term2 = (-3.0*b/2.0 + alpha/4.0) * ( 1.0 + (self.x/r_s)**2.0 ) ** (-3.0*b/2.0 + alpha/4.0 - 1.0) * 2.0/r_s * (self.x/r_s)
            d_term3 = -e/2.0 * ( 1.0 + (self.x/r_x)**3.0 ) ** (-e/6.0 - 1.0) * (self.x/r_x)**2.0 / r_x

            d_rho_inner = (d_term1 * term2 * term3) + (term1 * d_term2 * term3) + (term1 * term2 * d_term3)

        d_ftrans    = (1.0 + (self.x / r_t) ** beta) ** (-gamma / beta - 1.0) \
                      * (-gamma / r_t) * (self.x / r_t) ** (beta - 1.0)

        if self.outer_pro == 'MATTER':
            d_rho_outer = self.rho_mean * b_e * (-s_e) \
                          * (self.x / (5.0 * self.R200mean)) ** (-s_e - 1.0) \
                          * 1.0 / (5.0 * self.R200mean)
        elif self.outer_pro == 'BARYONIC':
            d_rho_outer = (self.Ob / self.Om) * self.rho_mean * b_e * (-s_e) \
                          * (self.x / (5.0 * self.R200mean)) ** (-s_e - 1.0) \
                          * 1.0 / (5.0 * self.R200mean)

        slope = d_rho_inner*ftrans + rho_inner*d_ftrans + d_rho_outer

        # derivatives of logs
        log_deriv = slope * self.x / profile
        if MODE == 'FIT':
            large_radii = np.where(self.x>.1)[0]
            #return log_deriv[large_radii] - self.profile_num_deriv[large_radii]
            return (log_deriv - self.profile_num_deriv)
        elif MODE == 'MODEL':
            return log_deriv

    def compute_splashback_radius(self):
        """
        Compute the splashback radius from the density derivative
        """

        #low_bound = .1 * self.R200mean

        #d2rho_dr2 = self.compute_numerical_deriv(self.x, self.fit_num_deriv, log=False)
        #print(d2rho_dr2)
        #print((d2rho_dr2 == np.max(d2rho_dr2)))
        #print((self.x > low_bound))
        #idx       = np.where((d2rho_dr2 == np.max(d2rho_dr2[(self.x > low_bound)])) )[0][0] + 1
        #Dint      = interp1d(d2rho_dr2[:idx], self.x[:idx], fill_value='extrapolate')
        #self.Rsp  = Dint(0.0)


        def deriv(r):

            # Original functions
            if self.inner_pro == 'EINASTO':
                rho_inner = self.rho_s * np.exp(-2.0 / self.alpha \
                                                * ((r / self.r_s) ** self.alpha - 1.0))
            elif self.inner_pro == 'BETA':
                rho_inner = self.rho_s / (1.0 + (r / self.r_s) ** 2.0) ** (1.5 * self.alpha)
            elif self.inner_pro == 'VIK06':
                rho_inner = self.rho_s * ((r / self.r_s) ** (-self.alpha * 0.5) \
                                    / (1.0 + (r / self.r_s) ** 2.0) \
                                          ** (3.0 * self.b / 2.0 - self.alpha / 4.0)) \
                    * (1.0 / ((1.0 + (r / self.r_x) ** 3.0) ** (self.e / 6.0)))

            ftrans    = (1.0 + (r / self.r_t) ** self.beta) ** (-self.gamma / self.beta)

            if self.outer_pro == 'MATTER':
                rho_outer = self.rho_mean * (self.b_e * (r / (5.0 * self.R200mean)) ** -self.s_e + 1.0)
            elif self.outer_pro == 'BARYONIC':
                rho_outer = (self.Ob / self.Om) * self.rho_mean * \
                            (self.b_e * (r / (5.0 * self.R200mean)) ** -self.s_e + 1.0)


            profile = rho_inner * ftrans + rho_outer

            # Gradients
            if self.inner_pro == 'EINASTO':
                d_rho_inner = -2.0 / self.r_s * (r / self.r_s) ** (self.alpha - 1.0) * rho_inner
            elif self.inner_pro == 'BETA':
                d_rho_inner = -3.0 * self.rho_s * self.alpha * r \
                              * ((1.0 + (r / self.r_s) ** 2.0) ** (-1.5 * self.alpha - 1.0)) \
                              / self.r_s ** 2.0
            elif self.inner_pro == 'VIK06':
                term1   = self.rho_s * (r/self.r_s) ** (-self.alpha * 0.5)
                term2   = ( 1.0 + (r/self.r_s)**2.0 ) ** (-3.0*self.b/2.0 + self.alpha/4.0)
                term3   = ( 1.0 + (r/self.r_x)**3.0 ) ** (-self.e/6.0)
                d_term1 = -self.alpha/2.0 * (r/self.r_s) ** (-self.alpha*0.5 - 1.0) * self.rho_s / self.r_s
                d_term2 = (-3.0*self.b/2.0 + self.alpha/4.0) * ( 1.0 + (r/self.r_s)**2.0 ) ** (-3.0*self.b/2.0 + self.alpha/4.0 - 1.0) * 2.0/self.r_s * (r/self.r_s)
                d_term3 = -self.e/2.0 * ( 1.0 + (r/self.r_x)**3.0 ) ** (-self.e/6.0 - 1.0) * (r/self.r_x)**2.0 / self.r_x

                d_rho_inner = (d_term1 * term2 * term3) + (term1 * d_term2 * term3) + (term1 * term2 * d_term3)


            d_ftrans    = (1.0 + (r / self.r_t) ** self.beta) ** (-self.gamma / self.beta - 1.0) \
                          * (-self.gamma / self.r_t) * (r / self.r_t) ** (self.beta - 1.0)

            if self.outer_pro == 'MATTER':
                d_rho_outer = self.rho_mean * self.b_e * (-self.s_e) \
                              * (r / (5.0 * self.R200mean)) ** (-self.s_e - 1.0) \
                              * 1.0 / (5.0 * self.R200mean)




            elif self.outer_pro == 'BARYONIC':
                 d_rho_outer = (self.Ob / self.Om) * self.rho_mean * self.b_e * (-self.s_e) \
                               * (r / (5.0 * self.R200mean)) ** (-self.s_e - 1.0) \
                               * 1.0 / (5.0 * self.R200mean)

            slope = d_rho_inner*ftrans + rho_inner*d_ftrans + d_rho_outer

            # Return complete
            return slope * r / profile

        lower_ind = 0
        upper_ind = -1
        while (lower_ind < (len(self.x)-1)) and (deriv(self.x[lower_ind]) < deriv(self.x[lower_ind+1])):
            lower_ind += 1
        if lower_ind == len(self.x)-1:
            lower_ind = 0
        while (upper_ind > (-len(self.x)+1)) and (deriv(self.x[upper_ind]) < deriv(self.x[upper_ind-1])):
            upper_ind -= 1
        if upper_ind == (-len(self.x)+1):
            upper_ind = -1
        try:
            Rguess = self.x[np.argmin(deriv(self.x[lower_ind:upper_ind]))+lower_ind]
        except:
            Rguess = self.R200mean

        bnds = [(self.x[lower_ind],self.x[upper_ind])]
        Rspl = minimize(deriv, Rguess, bounds=bnds)

        self.Rsp       = Rspl.x[0]
        self.Rsp_slope = deriv(self.Rsp)


        #del d2rho_dr2, idx, Dint
        return

    def density_fit_plot(self):
        """
        This function plots the profile and the analytic model fit
        """

        fig = plt.figure(figsize=(5,5))
        ax  = fig.add_subplot(111)
        pl1 = ax.loglog(self.x, self.y / self.rho_mean, 's', c='darkgrey', ms=4, \
                        label='Data', zorder=-4)
        if np.min(self.x) < 1.:
            ax.set_xlim(.009,5.5)
            ax.set_xlabel(r'$r/r_{200m}$', fontsize=18, family='serif')
        else:
            ax.set_xlim(1.0e1, 1.0e4)
            ax.set_xlabel(r'$r$  $[\mathrm{kpc}]$', fontsize=18, family='serif')
        ax.set_ylim(5.0e-5, 5.0e5)
        ax.set_yticks([1.0e-4, 1.0e-2, 1.0e0, 1.0e2, 1.0e4])
        ax.set_ylabel(r'$\rho\,/\,\rho_{\mathrm{mean}}$', fontsize=18, family='serif')
        ax.legend(loc=1, fontsize=10)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        fig.tight_layout()
        fig.savefig('images/{0}_density_fit.png'.format(self.tag), dpi=128)
        plt.close(fig)
        return

    def slope_comparison_plot(self):
        """
        This function plots the numerical derivative of the profile
        against the numerical derivative of the fit and the density
        derivative derived from the derivative of the analytic fit
        """

        fig = plt.figure(figsize=(5,5))
        ax  = fig.add_subplot(111)
        pl1 = ax.semilogx(self.x, self.profile_num_deriv, 's', c='darkgrey', ms=4, \
                          label='Num. deriv.', zorder=-4)
        pl2 = ax.semilogx(self.x, self.fit_num_deriv, ':', c='dodgerblue', lw=2, \
                          label='Fit num. deriv.', zorder=-2)
        pl3 = ax.semilogx(self.x, self.analytic_deriv, '-.', c='crimson', lw=2, \
                          label='Analytic deriv.', zorder=-1)
        pl4 = ax.semilogx([self.Rsp, self.Rsp], [-5.0, 2.0], '-', c='goldenrod', lw=2, \
                          label='Analytic Rsp')
        if hasattr(self, 'smooth_deriv'):
            pl5 = ax.semilogx(self.smooth_radius, self.smooth_deriv, '-.', c='black', lw=2, \
                              label='Smooth deriv.')
            pl6 = ax.semilogx([self.smooth_Rsp, self.smooth_Rsp], [-5.0,2.0],'-.', c='purple', lw=2, \
                              label='Smooth Rsp')
        pl7 = ax.semilogx([self.global_min_x, self.global_min_x], [-5.0, 2.0], '-.', c='green', lw=2, \
                          label='Global min.')
        if np.min(self.x) < 1.:
            ax.set_xlim(.009,5.5)
            ax.set_xlabel(r'$r/r_{200m}$', fontsize=18, family='serif')
        else:
            ax.set_xlim(1.0e1, 1.0e4)
            ax.set_xlabel(r'$r$  $[\mathrm{kpc}]$', fontsize=18, family='serif')
        ax.set_ylim(-5.0, 1.5)
        ax.set_ylabel(r'$d\log\rho\,/\,d\log r$', fontsize=18, family='serif')
        ax.legend(loc=2, fontsize=10)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        fig.tight_layout()
        fig.savefig('images/{0}_density_derivative.png'.format(self.tag), dpi=128)
        plt.close(fig)
        return




    def variance(self,ranged_var=False):

        if self.inner_pro in ['EINASTO', 'BETA']:
            p0 = [self.rho_s, self.r_s, self.alpha, self.r_t, self.beta, self.gamma, \
                  self.b_e, self.s_e]
        elif self.inner_pro in ['VIK06']:
            p0 = [self.rho_s, self.r_s, self.alpha, self.b, self.r_x, self.e, \
                  self.r_t, self.beta, self.gamma, self.b_e, self.s_e]
        if hasattr(self,'density_normalization'):
            fit_type = 'NORM'
        else:
            fit_type = 'MODEL'
        model = self.DK14_function(p0, self.inner_pro, self.outer_pro, fit_type)

        deriv = self.analytic_deriv

        model_var = np.sum( np.log10(model/self.y)**2 )
        deriv_var = np.sum( (deriv-self.profile_num_deriv)**2 )

        return_list = [model_var, deriv_var]

        if ranged_var:
            range_i = np.where( (self.x>(.8*self.R200mean)) & (self.x<(2.*self.R200mean)) )[0]
            ranged_var = np.sum( (deriv[range_i]-self.profile_num_deriv[range_i])**2 )
            return_list.append(ranged_var)

        return return_list



    def chi_squared(self, SLOPE=True):

        if self.inner_pro in ['EINASTO','BETA']:
            p0 = [self.rho_s, self.r_s, self.alpha, self.r_t, self.beta, self.gamma, \
                  self.b_e, self.s_e]
        elif self.inner_pro in ['VIK06']:
            p0 = [self.rho_s, self.r_s, self.alpha, self.b, self.r_x, self.e, \
                  self.r_t, self.beta, self.gamma, self.b_e, self.s_e]

        if self.selection is None:
            dof = 8
        else:
            dof = 5

        n_points_tot    = len(self.x)
        ranged_inds     = np.where( (self.x>(.8*self.R200mean)) & (self.x<(2.*self.R200mean)) )[0]
        n_points_ranged = len(ranged_inds)


        if SLOPE:
            model = self.analytic_deriv
            square_diff = (self.profile_num_deriv - model)**2
        else:
            model = self.DK14_function(p0, self.inner_pro, self.outer_pro, 'MODEL')
            square_diff = (self.y - model)**2

        chi_squared_tot = np.sum(square_diff)/(n_points_tot-dof)
        chi_squared_ranged = np.sum(square_diff[ranged_inds])/(n_points_ranged-dof)

        return chi_squared_tot, chi_squared_ranged




    def A14_Rsp(self, s):

        A = 38
        b = .57
        c = .02
        d = 0.2
        e = .52
        Ds = A * self.Om ** (-b - c*s) * np.exp(d*self.Om + e * s**(3./4.))

        Rint = interp1d(self.y / self.rho_mean, self.x)
        return Rint(Ds)

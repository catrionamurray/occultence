from ..imports import *


def lnprob(p,y,gp):
    # Trivial uniform prior.
    if np.any((-10 > p[1:]) + (p[1:] > 10)):
        return -np.inf

    # Update the kernel and compute the lnlikelihood.
    gp.set_parameter_vector(p)
    return gp.lnlikelihood(y, quiet=True)

def nllGP(p,y,t,e,gp):
    # if p[-1] < -1:
    #     return 1e25
    # print(p)
    gp.set_parameter_vector(p)
    # try:
    #     gp.compute(t, yerr=e)
    # except:
    #     return 1e25
    # c = -gp.log_likelihood(y)
    # return -gp.log_likelihood(y)
    ll = gp.log_likelihood(y, quiet=True)
    return -ll if np.isfinite(ll) else 1e25

def grad_nllGP(p,y,t,e,gp):
    gp.set_parameter_vector(p)
    # try:
    #     gp.compute(t, yerr=e)
    # except:
    #     return np.array([0.0]*len(p))
    # c=  -gp.grad_log_likelihood(y)
    return -gp.grad_log_likelihood(y,quiet=True)

def gp(x,
       y,
       yerr,
       x_original,
       rotation_period=None,
       rotation_amp=None,
       plot=False,
       figsize=(6,4),
       ):

    x_pred = np.linspace(np.min(x), np.max(x), 1000)

    # try:
    # jitter = george.kernels.ConstantKernel(log_constant=np.log(np.nanmedian(yerr)))
    amp = george.kernels.ConstantKernel(log_constant=np.log(np.std(y)))
    sqexp = george.kernels.ExpSquaredKernel(0.5, metric_bounds={'log_M_0_0':(np.log(0.01),np.log(1000))})

    # If the user has passed a rotation period to the GP then use quasi-periodic kernel, otherwise use a squared
    # exponential kernel.
    if rotation_period is not None:
        if rotation_amp is None:
            message = f""" A rotation period ({rotation_period}d), but no amplitude, has been passed to the GP, 
            therefore we will use 0.1 as the starting value"""
            cheerfully_suggest(message)
            rotation_amp = 0.1
        sinexp = george.kernels.ExpSine2Kernel(gamma=rotation_amp,log_period=np.log(rotation_period))
        gpkernel = amp * (sqexp*sinexp)
        kernel = "quasi-periodic"
        print(f"Fitting rotation with quasi-periodic GP, period = {rotation_period}d...")
    else:
        gpkernel = amp * sqexp
        kernel = "square_exponential"
        print(f"Fitting rotation with square-exponential GP...")

    gaussproc = george.GP(gpkernel,white_noise=np.log(np.nanmedian(yerr)),fit_white_noise=True)
    gaussproc.compute(x, yerr)

    print("Initial Parameter Vector: ", gaussproc.get_parameter_vector())
    if kernel == "quasi-periodic":
        p = gaussproc.get_parameter_vector()
        print('Initial Params: AMP = ', str(math.exp(p[1])), ', SQEXP= ', str(math.exp(p[2])),
              ', GAMMA= ', str(p[3]), ", Period = ", str(math.exp(p[4])), " days",', JITTER = ',
              str(math.exp(p[0])))
    else:
        p = gaussproc.get_parameter_vector()
        print('Initial Params: AMP = ', str(math.exp(p[1])), ', SQEXP= ', str(math.exp(p[2])),
              ', JITTER = ', str(math.exp(p[0])))

    print("Initial ln-likelihood: {0:.2f}".format(gaussproc.log_likelihood(y)))

    soln = minimize(nllGP, gaussproc.get_parameter_vector(), bounds=gaussproc.get_parameter_bounds(),
                    jac=grad_nllGP, method="L-BFGS-B", args=(y, x, yerr, gaussproc))
    p1 = soln.x

    print('Fitted GP HPs:', p1)
    if kernel == "quasi-periodic":
        print('Fitted Params: AMP = ', str(math.exp(p1[1])), ', SQEXP= ', str(math.exp(p1[2])), ', GAMMA= ',
              str(p1[3]), ", Period = ",
              str(math.exp(p1[4])), " days", ', JITTER = ', str(math.exp(p1[0])))
        hp = {'amp': math.exp(p1[1]), 'sqexp': math.exp(p1[2]), 'gamma': p1[3], 'period': math.exp(p1[4]),
              'jitter': math.exp(p1[0])}
    else:
        p = gaussproc.get_parameter_vector()
        print('Fitted Params: AMP = ', str(math.exp(p1[1])), ', SQEXP= ', str(math.exp(p1[2])), ', JITTER = ',
              str(math.exp(p1[0])))
        hp = {'amp': math.exp(p1[1]), 'sqexp': math.exp(p1[2]), 'jitter': math.exp(p1[0])}

    gaussproc.set_parameter_vector(p1)
    gaussproc.compute(x, yerr)
    jitter = math.exp(p[0])

    print("Final ln-likelihood: {0:.2f}".format(gaussproc.log_likelihood(y)))
    mu_zero, var = gaussproc.predict(y, x, return_var=True)
    mu = mu_zero + 1

    kernel_dict = {'name':kernel, 'hyperparameters':hp}

    if plot:
        plt.figure(figsize=figsize)
        plt.errorbar(x, y+1, yerr=yerr, fmt=".k", capsize=0,zorder=0, alpha=0.3)

        if x_pred is not None:
            pred_mu, pred_var = gaussproc.predict(y, x_pred, return_var=True)
            y_pred = pred_mu + 1
            yerr_pred = np.sqrt(pred_var)
        else:
            x_pred = x
            y_pred = mu
            yerr_pred = np.sqrt(var)

        plt.fill_between(x_pred, y_pred - yerr_pred, y_pred + yerr_pred,
                            color="orange", alpha=0.5,zorder=2)
        plt.plot(x_pred, y_pred, "orange", lw=1.5, alpha=0.8,zorder=2)
        plt.xlabel("Time")
        plt.ylabel("Relative Flux")
        plt.show()
        plt.close()

    # except Exception as e:
    #     print(e)
    #     return [],0, 0, {}
    og_mu, og_var = gaussproc.predict(y, x_original, return_var=True)

    return mu, var, jitter, og_mu, og_var, kernel_dict
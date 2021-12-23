# Function to compute the Information concatenate
# Based on: Li et al, 2012: https://academic.oup.com/biostatistics/article/13/3/523/248016
# See also: https://haroldpimentel.wordpress.com/2014/12/08/in-rna-seq-2-2-between-sample-normalization/
import numpy
import h5py
import scipy.signal


def compute_information_content(counts, control, alpha_low=0.05, alpha_high=0.5, smooth=True, smooth_stddev=250000, binsize=None):
    chroms = [c for c in counts]

    if smooth:
        assert binsize, 'For smoothing, a binsize needs to be provided.'
        W = scipy.signal.gaussian(numpy.ceil(smooth_stddev / binsize * 6 / 2) * 2 + 1, smooth_stddev / binsize)
        W /= W.sum()
        counts = {c: scipy.signal.fftconvolve(counts[c], W, mode='same') for c in chroms}
        control = {c: scipy.signal.fftconvolve(control[c], W, mode='same') for c in chroms}

    # generating matrix
    lin_control = numpy.concatenate([control[c] for c in chroms])
    lin_counts = numpy.concatenate([counts[c] for c in chroms])
    Xij = numpy.array([lin_counts, lin_control]).T
    Xij = Xij[(Xij > 0).sum(axis=1) > 1,]

    # setting initial condition
    S = numpy.ones(len(Xij), dtype=bool)
    Xi = Xij.sum(axis=1)
    sj = Xij[S].sum(axis=0) / Xi[S].sum(axis=0)

    sj_initial = sj.copy()

    # computing initial GOF --> power2 has been removed to force matching down
    GOFi = ( (Xij - numpy.outer(Xi, sj)) / (numpy.outer(Xi, sj)) ).sum(axis=1)

    assert len(GOFi) == len(S)

    # iteratively updating the bins that should be used
    _GOF_low, _GOF_high = numpy.percentile(GOFi, 100 * numpy.array([alpha_low, alpha_high]))

    for _ in range(20):
        S_update = (GOFi >= _GOF_low) & (GOFi < _GOF_high)

        if (S_update != S).sum() < 5:
            break

        S = S_update.copy()

        sj = Xij[S].sum(axis=0) / Xi[S].sum(axis=0)

        GOFi = ((Xij - numpy.outer(Xi, sj)) / (numpy.outer(Xi, sj))).sum(axis=1)

        _GOF_low, _GOF_high = numpy.percentile(GOFi, 100 * numpy.array([alpha_low, alpha_high]))


    sj_final = sj.copy()
    S_final = S.copy()
    GOFi_final = GOFi.copy()

    ratios = sj_final / sj_initial
    corr_factor = ratios[1] / ratios[0]
    return corr_factor

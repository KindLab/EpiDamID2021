# Custom scripts to compute region matrices and signal matrices to visualize enrichment heatmaps and average enrichment plots

import numpy
import matplotlib.pyplot
import h5py
import seaborn
import os

# function: calculate the bin coordinates for all regions
def computeMatrix_referencepoint(regions, binsize, before, after, refpoint='TSS', N='default'):
    '''
    The function computes a matrix of bin coordinates that can be used to efficiently extract
    sample counts/signal around a list of genomic coordinates. For each genomic locus in a region
    dataframe, the function calculates bin coordinates for sampling points around a reference point
    (e.g. TSS).

    Input:

    regions    a bedstyle dataframe (note the order of the columns in bed files)

    binsize    the size of the (equally sized) bins used in the sample data

    before     the length of the region (in bp) before the reference point

    after      the length of the region (in bp) after the reference point

    refpoint   ['TSS', 'TES', 'center'] The reference point on which to focus.

    N          The number of sampling points to be taken between 'before' and 'after'.
               Default is to take max(number of bins in the region, 50).

    Output:

    A dictionary containing per chromosome a numpy.array of dimensions (number of regions, N)
    containing the bin coordinates for each sampling point.

    Latest updated on 20200303

    '''

    # formatting regions dataframe
    regions = regions.reset_index(drop=True)
    assert regions.shape[1] >= 3, 'Invalid dataframe. Require at least [chrom, start, end].'
    if (regions.shape[1] < 6) and ('strand' not in regions.columns):
        regions = regions.iloc[:,:3]
        regions.columns = ['chr', 'start', 'end']
        regions.loc[:,'strand'] = '+'
    elif regions.shape[1] >= 6:
        regions = regions.iloc[:,:6]
        regions.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']

    # determining the steps relative to the reference point
    if N == 'default':
        region_length = before + after
        num_bins = int(numpy.ceil(region_length / binsize))
        N = num_bins if num_bins > 50 else 50
    steps = numpy.linspace(-before, after, N).astype(int)

    # compute the matrix
    region_matrix = {}

    for chrom, regions_chr in regions.groupby('chr'):
        chr_coordinates = list()

        if chrom.startswith('chr'):
            chrom = chrom.strip('chr')

        orig_indices = list() # required to restore original element order

        for strand, df in regions_chr.groupby('strand'):
            orig_indices.append(df.index.values)
            assert strand in ['+', '-', '.']

            if strand in ['+', '.']:
                if refpoint == 'TSS':
                    gene_coordinates = numpy.repeat(df['start'].values, N).reshape((len(df), N))
                elif refpoint == 'TES':
                    gene_coordinates = numpy.repeat(df['end'].values, N).reshape((len(df), N))
                elif refpoint == 'center':
                    mids = (df['start'] + df['end']).values/2
                    gene_coordinates = numpy.repeat(mids, N).reshape((len(df), N))
                gene_coordinates = ((gene_coordinates + steps[None,:]) // binsize).astype(int)
            elif strand == '-':
                if refpoint == 'TSS':
                    gene_coordinates = numpy.repeat(df['end'].values, N).reshape((len(df), N))
                elif refpoint == 'TES':
                    gene_coordinates = numpy.repeat(df['start'].values, N).reshape((len(df), N))
                elif refpoint == 'center':
                    mids = (df['start'] + df['end']).values/2
                    gene_coordinates = numpy.repeat(mids, N).reshape((len(df), N))
                gene_coordinates = ((gene_coordinates - steps[None,:]) // binsize).astype(int)

            chr_coordinates.append(gene_coordinates)

        orig_indices = numpy.concatenate(orig_indices)
        sort_ind = numpy.argsort(orig_indices)
        chr_coordinates = numpy.concatenate(chr_coordinates, axis=0)[sort_ind]

        region_matrix[chrom] = numpy.array(chr_coordinates)

    return region_matrix

# function: calculate the bin coordinates for all regions
def computeMatrix_scaleregions(regions, binsize, before, after, scale, N='default'):
    '''
    The function computes a matrix of bin coordinates that can be used to efficiently extract
    sample counts/signal around a list of genomic coordinates. For each genomic locus in a region
    dataframe, the function scales the region of interest (TSS-TES) and calculates bin coordinates
    for sampling points around and including the scaled region.

    Input:

    regions    a bedstyle dataframe (note the order of the columns in bed files)

    binsize    the size of the (equally sized) bins used in the sample data

    before     the length of the region (in bp) before the reference point

    after      the length of the region (in bp) after the reference point

    scale      the length to which all regions are scaled

    N          the total number of sampling steps to be taken in the 'before' and 'after' regions.
               The number of sampling steps for the scaled region will be proportional to the scale size.
               Default is to take max(number of bins in the region, 50).

    Output:

    A dictionary containing per chromosome a numpy.array of dimensions (number of regions on chr, N*(before + after + scale)/(before + after))
    containing the bin coordinates for each sampling point.

    Latest updated on 20200303

    '''

    # formatting regions dataframe
    regions = regions.reset_index(drop=True)
    assert regions.shape[1] >= 3, 'Invalid dataframe. Require at least [chrom, start, end].'
    if (regions.shape[1] < 6) and ('strand' not in regions.columns):
        regions = regions.iloc[:,:3]
        regions.columns = ['chr', 'start', 'end']
        regions.loc[:,'strand'] = '+'
    elif regions.shape[1] >= 6:
        regions = regions.iloc[:,:6]
        regions.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']

    # determining the steps relative to the reference point
    if N == 'default':
        num_bins = int(numpy.ceil(before / binsize))
        N_before = num_bins * 3
        num_bins = int(numpy.ceil(after / binsize))
        N_after = num_bins * 3
        if N_before + N_after < 150:
            f = numpy.ceil(150/(N_before + N_after))
            N_before = int(N_before * f)
            N_after = int(N_after * f)
        N = N_before + N_after
    else:
        N_before = int(numpy.round(N * before/(before + after), 0))
        N_after = int(numpy.round(N * after / (before + after), 0))
    N_scale = int(N * scale / (before + after))

    steps_before = numpy.linspace(-before, 0, N_before + 1).astype(int)[:-1] # the '0-coordinate' is not taken along, since this will be included in the scaled region
    steps_after = numpy.linspace(0, after, N_after + 1).astype(int)[1:] # the '0-coordinate' is not taken along, since this will be included in the scaled region

    # compute the matrix
    region_matrix = {}

    for chrom, regions_chr in regions.groupby("chr"):
        chr_coordinates = []

        if chrom.startswith('chr'):
            chrom = chrom.strip('chr')

        orig_indices = list()
        for strand, df in regions_chr.groupby("strand"):
            orig_indices.append(df.index.values)
            assert strand in ['+', '-', '.']

            # determine the steps within the scaled region
            if strand in ['+', '.']:
                starts = df['start'].values
                ends = df['end'].values
                b = numpy.repeat(starts, N_before).reshape((len(df), N_before)) + steps_before[None,:]
                a = numpy.repeat(ends, N_after).reshape((len(df), N_after)) + steps_after[None,:]
                s = numpy.repeat(starts, N_scale).reshape(len(starts), N_scale)
                reps = numpy.repeat(numpy.arange(N_scale), len(starts)).reshape(len(starts), N_scale, order='F')
                s = s + reps * ((ends - starts)/(N_scale-1))[:,None]

            elif strand == '-':
                starts = df['start'].values
                ends = df['end'].values
                b = numpy.repeat(ends, N_before).reshape((len(df), N_before)) - steps_before[None,:]
                a = numpy.repeat(starts, N_after).reshape((len(df), N_after)) - steps_after[None,:]
                s = numpy.repeat(ends, N_scale).reshape(len(ends), N_scale)
                reps = numpy.repeat(numpy.arange(N_scale), len(ends)).reshape(len(ends), N_scale, order='F')
                s = s + reps * ((starts - ends)/(N_scale-1))[:,None]

            gene_coordinates = numpy.concatenate([b,s,a], axis=1)
            gene_coordinates = (gene_coordinates // binsize).astype(int)
            chr_coordinates.append(gene_coordinates)


        orig_indices = numpy.concatenate(orig_indices)
        sort_ind = numpy.argsort(orig_indices)
        chr_coordinates = numpy.concatenate(chr_coordinates, axis=0)[sort_ind]
        region_matrix[chrom] = numpy.array(chr_coordinates)

    return region_matrix


# function: calculate all sample values for the
def computeSampleMatrix(sample_counts, region_matrix):
    '''
    The function calculates sample counts/signal for coordinates provided in a matrix computed by
    computeMatrix_referencepoint or computeMatrix_scaleregions.

    Input:

    sample_counts     A dictionary containing for each chromosome a list/numpy.array of counts/signal
                      observed in a bin.

    region_matrix     A dictionary containing for each chromosome an numpy.array of bin coordinates.
                      Output of computeMatrix_referencepoint or computeMatrix_scaleregions.

    Output:

    A matrix of shape (number of regions, number of steps + 1) containing the signal observed in the bins.

    Last updated on 20200303

    '''
    sample_matrix = list()

    for chrom, chr_counts in sample_counts.items():
        if chrom not in region_matrix:
            continue
        chr_coords = region_matrix[chrom]
        chr_coords[chr_coords < 0] = 0
        chr_coords[chr_coords >= len(chr_counts)] = len(chr_counts) - 1

        sample_matrix.append(chr_counts[chr_coords])

    sample_matrix = numpy.concatenate(sample_matrix, axis=0)

    return sample_matrix


def plotHeatmap(sample_data, before, after, regiontype='TSS', scale_length=None, plotlen=10, cmaprange='default', yrange='default'):
    '''
    The function plots heatmaps and summary plots based on matrices produced by
    computeMatrix_scaleregions/computeMatrix_referencepoint and computeSampleMatrix.

    Input:

    sample_data       A matrix outputted by computeSampleMatrix or a dictionary containing for several
                      samples such a matrix.

    before            The length of the region upstream of the reference point / scaled region

    after             The length of the region downstream of the reference point / scaled region

    regiontype        ['TSS', 'TES', 'center', 'scale'] The type of regions plotted.

    scale_length      if regiontype == 'scale', the length to which all regions are scaled is required.

    plotlen           The length of the produced figure. default = 10.

    cmaprange         The range of the values used for the heatmap color bar.

    yrange            the range of the values used for the y-axis in the summary plot.

    Output:

    A figure showing the summary plot and heatmaps of the samples. Use matplotlib.pyplot.show() to visualize.

    Last updated on 180423

    '''

    if regiontype == 'scale':
        assert scale_length is not None, 'Region type \'scale\' requires a scale length'
        assert type(scale_length) == int, 'Scale length has to be an integer'
        total_len = before + after + scale_length

    if type(sample_data) == dict:
        to_plot = sample_data
    else:
        to_plot = {}
        to_plot['sample'] = sample_data

    fig = matplotlib.pyplot.figure(figsize=(2*len(to_plot.keys()), plotlen))
    row_elements = int(plotlen/2)
    column_elements = len(to_plot.keys())

    counter = 0

    for sample in sorted(list(to_plot.keys())):
        matrix = to_plot[sample]
        num_steps = matrix.shape[1]

        ### summary plot ###
        ax0 = matplotlib.pyplot.subplot2grid((row_elements, column_elements), (0, counter))
        x = numpy.arange(0, num_steps)
        y = matrix.mean(axis=0)
        ax0.plot(x, y, lw=2)

        if regiontype in ['TSS', 'TES', 'center']:
            ax0.axvline(x=(before/(before+after))*x[-1], c='black', alpha=0.5)
        elif regiontype == 'scale':
            s = before/total_len*x[-1]
            e = (before+scale_length)/total_len*x[-1]
            ax0.axvline(x=s, c='black', alpha=0.5)
            ax0.axvline(x=e, c='black', alpha=0.5)

        #formatting
        if yrange != 'default':
            ax0.set_ylim((yrange[0],yrange[1]))
        ax0.set_title(sample)

        ax0.set_facecolor('white')
        for sp in ax0.spines:
            ax0.spines[sp].set_color('black')
            ax0.spines[sp].set_linewidth(1)

        ### heatmap ###
        ax1 = matplotlib.pyplot.subplot2grid((row_elements, column_elements), (1, counter), rowspan=(row_elements - 1))
        sort_indices = numpy.argsort(matrix.sum(axis=1))

        if cmaprange == 'default':
            p1 = ax1.pcolorfast(matrix[sort_indices,:], cmap='Blues')
        else:
            p1 = ax1.pcolorfast(matrix[sort_indices,:], cmap='Blues', vmin=cmaprange[0], vmax=cmaprange[1])
        matplotlib.pyplot.colorbar(p1, ax=ax1, orientation='horizontal')

        #formatting
        ax1.set_yticks([0])

        ### formatting x labels
        offset = int(0.05*num_steps)
        for ax in [ax0, ax1]:
            ax.grid(False)
            if regiontype in ['TSS', 'TES', 'center']:
                ax.set_xticks([0 + offset, num_steps/2, num_steps - offset])
                ax.set_xticklabels(["-{}kb".format(before/1000), regiontype, "+{}kb".format(after/1000)], rotation='horizontal')
            elif regiontype == 'scale':
                ax.set_xticks([0 + offset, int(before*num_steps/total_len), int((before+scale_length)*num_steps/total_len), num_steps - offset])
                ax.set_xticklabels(["-{}kb".format(before/1000), "TSS", "TES", "+{}kb".format(after/1000)], rotation='horizontal')

        counter += 1

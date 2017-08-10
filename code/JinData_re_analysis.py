from matplotlib import pyplot as plt
from csv import reader as csv_reader
import numpy as np
import os
import matplotlib
from chiffatools.linalg_routines import gini_coeff, rm_nans
from csv import reader, writer
import numpy as np
from chiffatools.linalg_routines import gini_coeff
from matplotlib import pyplot as plt
from scipy.stats import expon

# source_folder = 'L:\\Users\\andrei\\Jin_2010_re_analysis'
# values = 'speeds_v.csv'
# errors = 'speeds_err.csv'

source_folder = 'JinData'
values = 'lastvals_v.csv'

# source_folder = 'JinDataGBO'
# values = 'meanNormalizedFinalAverageSpotIntensity.csv'

errors = 'speeds_err.csv'
outfile = 'JinData_out_haploid-gen.csv'
transpose = True


def csv2numpy(source, c_header=True, r_header=True):

    def correct_line(_row):
        return [float(item) if item not in ['inf', '', ' ', 'NA'] else np.inf for item in _row]

    with open(source, 'r') as source_file:
        reader = csv_reader(source_file)

        if c_header:
            c_headers = reader.next()
        else:
            c_headers = []

        r_headers = []
        data_container = []

        for row in reader:
            if r_header:
                r_headers.append(row[0])
                row = row[1:]
            data_container.append(correct_line(row))

        return np.array(data_container), c_headers[1:], r_headers

vals, c_headers, r_headers = csv2numpy(os.path.join(source_folder, values))
# errs, _, _ = csv2numpy(os.path.join(source_folder, errors))

synthetic = vals
# synthetic = 1/vals  # inversion from doubling time to growth speed

# print vals
# print synthetic


def clear_empty_rows_and_cols(_synthetic, c_headers, r_headers):
    _synthetic = _synthetic / np.max(np.ma.masked_invalid(_synthetic))
    # exclusion = ['low glucose', 'glycerol']

    gini_indexes = np.apply_along_axis(gini_coeff, 0, _synthetic)
    keep_cols = gini_indexes < 0.7

    # for i, name in enumerate(c_headers):
    #     if name in exclusion:
    #         keep_cols[i] = False

    _synthetic = _synthetic[:, keep_cols]
    c_headers = np.array(c_headers)[keep_cols]

    gini_indexes = np.apply_along_axis(gini_coeff, 1, _synthetic)
    keep_rows = gini_indexes < 0.7
    _synthetic = _synthetic[keep_rows, :]
    r_headers = np.array(r_headers)[keep_rows]

    _synthetic[_synthetic == 0] = np.nan

    # drop the first 3:
    _synthetic = _synthetic[3:, :]
    r_headers = np.array(r_headers)[3:]

    return _synthetic, c_headers, r_headers


synthetic, c_headers, r_headers = clear_empty_rows_and_cols(synthetic, c_headers, r_headers)

# adding noise
scale = np.min(synthetic[synthetic > 0])/expon.interval(0.95)[1]
synthetic_w_noise = synthetic.copy()
synthetic_w_noise[synthetic_w_noise == 0] = np.random.exponential(scale, size=synthetic.shape)[
    synthetic_w_noise == 0]

synthetic[np.isinf(synthetic)] = np.nan

# removing the haploids:
def analytic_color(matrix_of_interest):
    x_gini_indexes = np.apply_along_axis(gini_coeff, 0, matrix_of_interest)
    x_means = np.apply_along_axis(np.nanmean, 0, matrix_of_interest)
    x_pad = np.zeros_like(x_means)

    y_gini_indexes = np.apply_along_axis(gini_coeff, 1, matrix_of_interest)
    y_means = np.apply_along_axis(np.mean, 1, matrix_of_interest)
    y_pad = np.zeros_like(y_means)

    ret_matrix = np.vstack((np.array([x_gini_indexes, x_means, x_pad]),
                                     matrix_of_interest))
    ret_matrix = np.hstack((
        np.vstack((np.zeros((3, 3)), np.array([y_gini_indexes, y_means, y_pad]).T)),
        ret_matrix
                           ))

    return ret_matrix


def show_cleared_matrix():
    a_synthetic = analytic_color(synthetic)
    a_synthetic_w_noise = analytic_color(synthetic_w_noise)

    plt.subplot(121)
    plt.title('Original')
    plt.imshow(a_synthetic, interpolation='nearest')
    plt.xlabel('stresses')
    plt.ylabel('aneuploids')
    plt.colorbar()

    plt.subplot(122)
    plt.title('pre-processed')
    plt.imshow(a_synthetic_w_noise, interpolation='nearest')
    plt.xlabel('stresses')
    plt.ylabel('aneuploids')
    plt.colorbar()

    plt.show()

show_cleared_matrix()


def plot_raw_values():

    cmap = matplotlib.cm.get_cmap('Dark2')
    l = synthetic.shape[0]

    # argsort = np.argsort(synthetic[-2])[-1:0:-1]
    argsort = np.argsort(np.nanmean(synthetic, axis=0))[::-1]
    for i, (row, r_name) in enumerate(zip(np.vsplit(synthetic, synthetic.shape[0]), r_headers)):
        if r_name in ['U1', 'controlHaploid']:
            plt.plot(row[0, argsort], 'k')
            plt.plot(row[0, argsort], 'ko', label=r_name)
        else:
            plt.plot(row[0, argsort], 'o', label=r_name, color=cmap(float(i)/l))

    plt.xticks(np.arange(len(c_headers)), c_headers[argsort], rotation="vertical")
    plt.legend(ncol=2)
    plt.show()

    gini_indexes = np.apply_along_axis(gini_coeff, 1, synthetic)
    means = np.apply_along_axis(np.nanmean, 1, synthetic)

    print gini_indexes
    print means

    plt.title('Gini - mean fitness correlation')
    plt.plot(gini_indexes, means, 'ko', label="aneuploids")
    # print r_headers
    # print r_headers[38]
    # plt.plot(gini_indexes[38], means[38], 'ro', label="haploids")
    # print r_headers[[39, 43]]
    plt.plot(gini_indexes[[39, 43]], means[[39, 43]], 'ro', label="haploids")
    # print r_headers[[40, 42]]
    # plt.plot(gini_indexes[[40, 42]], means[[40, 42]], 'bo', label="diploids")
    # print r_headers[[41, 44]]
    # plt.plot(gini_indexes[[41, 44]], means[[41, 44]], 'go', label="triploids")
    plt.legend()
    plt.xlabel('Gini index')
    plt.ylabel('fitness')
    plt.show()


plot_raw_values()

# switching to stresses as aneuploids
if transpose:
    synthetic = synthetic.T



# print r_headers
# print r_headers[[39, 43]]
# print c_headers

gini_indexes = np.apply_along_axis(gini_coeff, 0, synthetic)[:-6]
means = np.apply_along_axis(np.nanmean, 0, synthetic)[:-6]

gini_ranks = np.arange(synthetic.shape[1])[np.argsort(gini_indexes)]
mean_ranks = np.arange(synthetic.shape[1])[np.argsort(means)[::-1]]

# combined_argsort = np.argsort(np.min(np.array([gini_ranks, mean_ranks]), axis=0))[::-1]
combined_argsort = np.argsort(gini_indexes)

# print gini_indexes[combined_argsort]
# print means[combined_argsort]

# TODO: exclude the pre-set amount of elements

generalist_inclusion = 3
synthetic_generalist = np.mean(synthetic[:, combined_argsort[0:generalist_inclusion]], axis=1)
running_lowest_gini = gini_coeff(synthetic_generalist)

while gini_coeff(synthetic_generalist) <= running_lowest_gini + 0.0005:
    running_lowest_gini = gini_coeff(synthetic_generalist)
    generalist_inclusion += 1
    synthetic_generalist = np.nanmean(synthetic[:, combined_argsort[0:generalist_inclusion]],
                                      axis=1)
    print generalist_inclusion, running_lowest_gini, gini_coeff(synthetic_generalist)

generalist_inclusion -= 1
selected = synthetic[:, combined_argsort[0:generalist_inclusion]]
print 'selected bases for synth generalist:', r_headers[combined_argsort[0:generalist_inclusion]]
# selected = synthetic[:, [42, 46]]
# TODO: there is a need to adjust the correction factor in a unique manner, so that it isnt' too
# high
synthetic_generalist = np.nanmean(synthetic[:, combined_argsort[0:generalist_inclusion]],
                                  axis=1) + 0.05
old_generalist = np.nanmean(synthetic[:, [39, 43]], axis=1)
# old_generalist = np.nanmean(synthetic[:, [38]], axis=1)
synthetic_generalist = old_generalist
running_lowest_gini = gini_coeff(synthetic_generalist)

# override for the actual generalist:

plt.title('raw plots')
plt.plot(np.log10(synthetic_generalist), 'ro', label="synthetic generalist")
plt.plot(np.log10(synthetic), 'k.')
plt.plot(np.log10(selected), 'b.', label="bases for synthetic generalist")
plt.legend(numpoints=1)
plt.xlabel('stresses')
plt.ylabel('fitness')
plt.xticks(np.arange(len(c_headers)), c_headers, rotation="vertical")
plt.show()

plt.title('raw plots')
plt.plot(np.log10(synthetic_generalist), 'ro', label="synthetic generalist")
plt.plot(np.log10(old_generalist), 'go', label="haploid")
plt.legend(numpoints=1)
plt.xlabel('odorants')
plt.ylabel('fitness')
plt.xticks(np.arange(len(c_headers)), c_headers, rotation="vertical")
plt.show()

gini_indexes = np.apply_along_axis(gini_coeff, 0, synthetic)
means = np.apply_along_axis(np.nanmean, 0, synthetic)

print gini_indexes
print means

plt.title('')
plt.plot(gini_coeff(synthetic_generalist), np.nanmean(synthetic_generalist), 'ro',
         label="synthetic generalist")
plt.plot(gini_indexes, means, 'k.')
plt.plot(np.apply_along_axis(gini_coeff, 0, selected),
         np.apply_along_axis(np.nanmean, 0, selected),
         'b.', label="bases for synthetic generalist")
# plt.plot(np.apply_along_axis(gini_coeff, 0, synthetic[:, [42, 46]]),
#          np.apply_along_axis(np.mean, 0, synthetic[:, [42, 46]]),
#          'g.', label="haploid")
plt.legend()
plt.xlabel('Gini index')
plt.ylabel('fitness')
plt.show()
# plt.close()

_norm = np.apply_along_axis(lambda x: x/synthetic_generalist, 0, synthetic)
print _norm

log2_norm = np.log2(_norm)
mean = np.apply_along_axis(lambda x: np.nanmean(x), 1, log2_norm)
std = np.apply_along_axis(lambda x: np.nanstd(x), 1, log2_norm)
mean_argsort = np.argsort(mean)[::-1]

print log2_norm

log2_norm = log2_norm[mean_argsort, :]
mean = np.apply_along_axis(lambda x: np.nanmean(x), 1, log2_norm)
std = np.apply_along_axis(lambda x: np.nanstd(x), 1, log2_norm)
c_headers = c_headers[mean_argsort]

plt.title('log2-normalized to synthetic generalist')
plt.plot(log2_norm, 'k.')
plt.plot(mean, 'r')
plt.plot(mean + std, 'g')
plt.plot(mean - std, 'g')
plt.xticks(np.arange(len(c_headers)), c_headers, rotation="vertical")
plt.show()

plt.title('mean-log correlation')
plt.plot(mean, std, 'k.')
plt.xlabel('n_log2 - Mean')
plt.ylabel('n_log2 - Standard deviation')
# plt.set_autoscalex_on(False)
# plt.set_autoscaley_on(False)
# plt.set_ylim([0, 1])
# plt.set_xlim([0, 1])
plt.show()

with open(outfile, 'wb') as destination:
    csv_writer = writer(destination)
    for mean, std in zip(mean.tolist(), std.tolist()):
        csv_writer.writerow([mean, std])

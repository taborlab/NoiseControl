# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 10:13:49 2017

@author: Evan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from copy import deepcopy
from openpyxl import load_workbook

class FlowSampleDataFrame():

    def __init__(self, data_file, sample_ids, channel):
        data_file = pd.ExcelFile(data_file)

        # Load and organize sample df
        sample_df = pd.read_excel(data_file, sheet_name='Samples') # sheetname -> sheet_name, 191226 kpg
        sample_df = sample_df[sample_df['ID'].isin(sample_ids)]

        # Load and organize bins and events arrays
        hist_df = pd.read_excel(data_file, sheet_name='Histograms') # sheetname -> sheet_name, 191226 kpg
        hist_df = hist_df[hist_df['Channel'] == channel]
        hist_df = hist_df[hist_df['Sample ID'].isin(sample_ids)]
        bin_labels = ['Bin {}'.format(i+1) for i in range(1024)]

        bins_df = hist_df[hist_df['Unnamed: 2'] == 'Bin Centers (MEF)']
        bins = bins_df[bin_labels].values

        hists_df = hist_df[hist_df['Unnamed: 2'] == 'Counts']
        hists = hists_df[bin_labels].values
        
        events = [hist_to_eventlist(bins[i], hists[i]) \
                    for i in range(len(bins))]

        # Set up a dictionary to map Sample ID to bins/events array indices
        # This is probably a more-robust-than-necessary approach.
        ID_to_index = dict(zip(hists_df['Sample ID'].values,
                          range(len(hists_df['Sample ID'].values))))

        index_df = pd.DataFrame.from_dict(ID_to_index, orient='index')
        index_df.rename(columns={0:'Data index'}, inplace=True)
        index_df.index.names = ['ID']

        # Merge sample_df and index_df by sample 'ID'
        sample_df.set_index('ID', inplace=True, drop=True)
        sample_df = pd.concat([sample_df, index_df], axis=1).reset_index()
        sample_df.rename(columns={'index':'ID'}, inplace=True)

        self.sample_df = sample_df
        self.bins = bins
        self.events = events
        self.channel = channel

def import_hists(fc_data_path, sample_ids, channel):
    return FlowSampleDataFrame(fc_data_path, sample_ids, channel)

def hist_to_eventlist(bin_centers, bin_counts):
    x = [[bin_centers[i]]*int(bin_counts[i]) for i in range(len(bin_centers)) if bin_counts[i] != 0]
    x = [item for sublist in x for item in sublist]
    return x

def kde_sklearn(x, x_grid, bandwidth=0.2, rtol=1e-4, **kwargs):
    """Kernel Density Estimation with Scikit-Learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return [np.exp(log_pdf), kde_skl]

def pdf_thresh(eventlist, thresh=0.01, bandwidth=0.05):
    log10_s1 = np.log10(eventlist)
    bins = np.linspace(0, 8, 400 + 1)
    s1_log_pdf, s1_log_kde = kde_sklearn(log10_s1, x_grid=bins, bandwidth=bandwidth, rtol=1e-4)

    prob_thresh = thresh * np.max(s1_log_pdf)
    median = np.median(log10_s1)

    bin_ints = np.linspace(0,400,401)
    idx = np.where(bin_ints == int(round(median*50)))[0][0]

    func = s1_log_pdf - prob_thresh

    # find low cutoff
    low_idx = idx
    while func[low_idx] > 0:
        low_idx -= 1

    # find high cutoff
    high_idx = idx
    while func[high_idx] > 0:
        high_idx += 1

    log10_s1_corr = log10_s1[np.where((log10_s1 > bins[low_idx]) & (log10_s1 < bins[high_idx]))]
    s1_corr = 10**log10_s1_corr

    return s1_corr, [bins[low_idx], bins[high_idx]]

def trim(fsdf, thresh=0.005, bw=0.05):
    fsdf_trimmed = deepcopy(fsdf)
    for idx, row in fsdf.sample_df.iterrows():
        data_index = int(row['Data index'])
        data_raw = fsdf.events[data_index]
        data_trimmed, trim_vals = pdf_thresh(data_raw, thresh, bw)

        fsdf_trimmed.events[data_index] = data_trimmed

        m = np.mean(data_trimmed)
        cv = np.std(data_trimmed)/m

        fsdf_trimmed.sample_df.loc[idx, '{} Trim Threshold'.format(fsdf.channel)] = thresh
        fsdf_trimmed.sample_df.loc[idx, '{} Trim Bandwidth'.format(fsdf.channel)] = bw
        fsdf_trimmed.sample_df.loc[idx, '{} Mean (Trimmed)'.format(fsdf.channel)] = m
        fsdf_trimmed.sample_df.loc[idx, '{} CV (Trimmed)'.format(fsdf.channel)] = cv
        fsdf_trimmed.sample_df.loc[idx, '{} Trim Low Result'.format(fsdf.channel)] = 10**trim_vals[0]
        fsdf_trimmed.sample_df.loc[idx, '{} Trim High Result'.format(fsdf.channel)] = 10**trim_vals[1]
        fsdf_trimmed.sample_df.loc[idx, '{} Number of Events (Trimmed)'.format(fsdf.channel)] = \
                                                    int(len(data_trimmed))
    return fsdf_trimmed

def plot(trimmed_fsdf, plot_folder, prefix=None, expt_id=None):
    # Make plots
    for idx, row in trimmed_fsdf.sample_df.iterrows():
        samp = row['ID']
        data_index = int(row['Data index'])
        m1 = row['{} Mean'.format(trimmed_fsdf.channel)]
        cv1 = row['{} CV'.format(trimmed_fsdf.channel)]
        m2 = row['{} Mean (Trimmed)'.format(trimmed_fsdf.channel)]
        cv2 = row['{} CV (Trimmed)'.format(trimmed_fsdf.channel)]
        len1 = int(row['Number of Events'])
        len2 = int(row['{} Number of Events (Trimmed)'.format(trimmed_fsdf.channel)])

        x = trimmed_fsdf.events[data_index]
        b = trimmed_fsdf.bins[data_index]
        trim_vals = [row['{} Trim Low Result'.format(trimmed_fsdf.channel)], row['{} Trim High Result'.format(trimmed_fsdf.channel)]]

        plt.clf()
        h = plt.hist(x, b[::10])
        plt.xscale('log')
        plt.plot([trim_vals[0]]*2, [0, np.max(h[0])], 'k-')
        plt.plot([trim_vals[1]]*2, [0, np.max(h[0])], 'k-')
        plt.plot([np.min(x)]*2, [0, np.max(h[0])], 'r-')
        plt.plot([np.max(x)]*2, [0, np.max(h[0])], 'r-')
        if expt_id:
            plt.title('Trimming - {} - {} \n #: {:d} ({:d}), Mean: {:.1f} ({:.1f}), CV: {:.3f} ({:.3f})'.format(
                       expt_id, samp, len2, len1, m2, m1, cv2, cv1))
        else:
            plt.title('Trimming - {} \n #: {:d} ({:d}), Mean: {:.1f} ({:.1f}), CV: {:.3f} ({:.3f})'.format(
                       samp, len2, len1, m2, m1, cv2, cv1))
        plt.xlim([1,1e6])
        plt.xlabel('Fluorescence ({})'.format(trimmed_fsdf.channel))
        plt.ylabel('Counts (#)')
        plt.subplots_adjust(top=0.82, right=0.97)
        if prefix:
            plt.savefig('{}/{}_{}.png'.format(plot_folder, prefix, samp))
        else:
            plt.savefig('{}/{}.png'.format(plot_folder, samp))

def af_correct_function(mean_tot, cv_tot, mean_af, cv_af):
    mean_fp = mean_tot - mean_af
    if mean_fp != 0:
        cv_fp = np.sqrt(cv_tot**2 * mean_tot**2 - cv_af**2 * mean_af**2) \
                                                        / (mean_tot - mean_af)
    else:
        cv_fp = np.nan
        
    return mean_fp, cv_fp

def af_correct(trimmed_fsdf, trimmed_fsdf_white):
    afcorr_fsdf = deepcopy(trimmed_fsdf)
    
    if trimmed_fsdf.channel == 'FL1':
        fluorophore = 'GFP'
    elif trimmed_fsdf.channel == 'FL3':
        fluorophore = 'RFP'
    else:
        fluorophore = 'FP'
    
    # Add af-corrected stats to trimmed_data
    mean_af = np.mean(trimmed_fsdf_white.sample_df['{} Mean'.format(trimmed_fsdf.channel)].values)
    cv_af = np.mean(trimmed_fsdf_white.sample_df['{} CV'.format(trimmed_fsdf.channel)].values)

    mean_af_trimmed = np.mean(trimmed_fsdf_white.sample_df['{} Mean (Trimmed)'.format(trimmed_fsdf.channel)].values)
    cv_af_trimmed = np.mean(trimmed_fsdf_white.sample_df['{} CV (Trimmed)'.format(trimmed_fsdf.channel)].values)

    for idx, row in trimmed_fsdf.sample_df.iterrows():
        mean_tot = row['{} Mean'.format(trimmed_fsdf.channel)]
        cv_tot = row['{} CV'.format(trimmed_fsdf.channel)]
        mean_fp, cv_fp = af_correct_function(
                    mean_tot, cv_tot, mean_af, cv_af)

        mean_tot_trimmed = row['{} Mean (Trimmed)'.format(trimmed_fsdf.channel)]
        cv_tot_trimmed = row['{} CV (Trimmed)'.format(trimmed_fsdf.channel)]
        mean_fp_trimmed, cv_fp_trimmed = af_correct_function(
                    mean_tot_trimmed, cv_tot_trimmed, mean_af_trimmed, cv_af_trimmed)

        afcorr_fsdf.sample_df.loc[idx, '{} Mean'.format(fluorophore)] = mean_fp
        afcorr_fsdf.sample_df.loc[idx, '{} CV'.format(fluorophore)] = cv_fp
        afcorr_fsdf.sample_df.loc[idx, '{} Mean (Trimmed)'.format(fluorophore)] = mean_fp_trimmed
        afcorr_fsdf.sample_df.loc[idx, '{} CV (Trimmed)'.format(fluorophore)] = cv_fp_trimmed

    return afcorr_fsdf

def save(trimmed_data, fc_data_path):
    # Save trimmed_data into new Excel Sheet in file at fc_data_path
    book = load_workbook(fc_data_path)
    writer = pd.ExcelWriter(fc_data_path, engine='openpyxl')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

    trimmed_data.sample_df.to_excel(writer, sheet_name='Trimmed')

    writer.save()

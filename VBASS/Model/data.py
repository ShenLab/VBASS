
import numpy as np
import torch
import scipy.sparse
import pandas as pd


class DataSet(torch.utils.data.IterableDataset):
    """
        Dataset object for loading TADA and gene expression input to autoencoder

        Parameters
        ----------
        dnv_table: pd.DataFrame
            An object of Genes with their mutation rates and de novo variants, column names
            should be named as cls1, cls2, etc.
        gene_list: pd.DataFrame
            A gene list of interesting genes, index should match the index names of DNV_table
            Other columns will be added to final output
            Optional:
            Can have a column named 'label' for semi-supervised training, otherwise unsupervised
        gene_scores: pd.DataFrame
            A data frame of gene scores to fed into autoencoder, index should match DNV_table
        shuffle: boolean
            Whether the data should be shuffled
        device: str
            cpu or cuda

        Parameters to be add
        ----------
        hypothesis_table: a table to indicate hypothesises. For now only two hypothesis
        """

    def __init__(
            self,
            # dnv_table: pd.DataFrame,
            # gene_list: pd.DataFrame,
            # gene_scores: pd.DataFrame,
            # cell_type_number: int,
            # shuffle=True,
            x,
            x_dim,
            label,
            x_var,
            train_flag: bool = True,
            warmup_flag: bool = True
    ):
        # gene_list = gene_list.iloc[np.isin(gene_list.index, dnv_table.index), :]
        # gene_list = gene_list.iloc[np.isin(gene_list.index, gene_scores.index), :]
        # dnv_table = dnv_table.loc[gene_list.index, :]
        # gene_scores = gene_scores.loc[gene_list.index, :]
        # index = np.arange(gene_list.shape[0])
        # if shuffle:
        #     np.random.shuffle(index)
        # gene_list = gene_list.iloc[index, :]
        # dnv_table = dnv_table.iloc[index, :]
        # gene_scores = gene_scores.iloc[index, :]

        # self.x_dim = [gene_scores.shape[1], dnv_table]

        # increase labeled data points
        if train_flag:
            x_labeled_index = np.where(~np.isnan(label.iloc[:, 0].values))[0]
            x_unlabeled_index = np.where(np.isnan(label.iloc[:, 0].values))[0]
            repeat_times = np.ceil(x_unlabeled_index.shape[0]/x_labeled_index.shape[0])
            x_labeled_index = np.repeat(x_labeled_index, repeats=repeat_times)
            x_index = np.concatenate((x_labeled_index, x_unlabeled_index))
            x = x.iloc[x_index, :]
            x_var = x_var.iloc[x_index, :]
            label = label.iloc[x_index, :]
        # shuffle x
        index = np.arange(x.shape[0])
        if train_flag:
            np.random.seed(0)
            np.random.shuffle(index)
        x = x.iloc[index, :]
        x_var = x_var.iloc[index, :]
        label = label.iloc[index, :]

        self.x = x
        self.x_dim = x_dim
        self.x_var = x_var
        self.label = label
        self.train_flag = train_flag

    def __iter__(self):
        """
        Iterates over all genes in the gene_list: self.gene_list in the data,
        yielding:
            x: node feature, in this case gene counts, with shape [num_nodes, 1]
            edge_index: single cell Graph connectivity in COO format with shape [2, num_edges]
            edge_attr: single cell Graph connectivity in COO format with shape [num_edges, 1]
            y: Graph classification
        """
        for i in range(self.x.shape[0]):
            yield (self.x.iloc[i, :].values.astype(np.float32),
                   self.x_var.iloc[i, :].values.astype(np.float32),
                   self.label.iloc[i, :].values.astype(np.float32)
                   if self.train_flag else np.array([np.nan]))

    def __len__(self):
        return len(self.x.shape[0])


class DataSet_v3(torch.utils.data.IterableDataset):
    """
        Dataset object for loading TADA and gene expression input to autoencoder

        Parameters
        ----------
        dnv_table: pd.DataFrame
            An object of Genes with their mutation rates and de novo variants, column names
            should be named as cls1, cls2, etc.
        gene_list: pd.DataFrame
            A gene list of interesting genes, index should match the index names of DNV_table
            Other columns will be added to final output
            Optional:
            Can have a column named 'label' for semi-supervised training, otherwise unsupervised
        gene_scores: pd.DataFrame
            A data frame of gene scores to fed into autoencoder, index should match DNV_table
        shuffle: boolean
            Whether the data should be shuffled
        device: str
            cpu or cuda

        Parameters to be add
        ----------
        hypothesis_table: a table to indicate hypothesises. For now only two hypothesis
        """

    def __init__(
            self,
            # dnv_table: pd.DataFrame,
            # gene_list: pd.DataFrame,
            # gene_scores: pd.DataFrame,
            # cell_type_number: int,
            # shuffle=True,
            x_in,
            x_out,
            x_dim,
            label,
            x_var,
            train_flag: bool = True,
            warmup_flag: bool = True
    ):
        # gene_list = gene_list.iloc[np.isin(gene_list.index, dnv_table.index), :]
        # gene_list = gene_list.iloc[np.isin(gene_list.index, gene_scores.index), :]
        # dnv_table = dnv_table.loc[gene_list.index, :]
        # gene_scores = gene_scores.loc[gene_list.index, :]
        # index = np.arange(gene_list.shape[0])
        # if shuffle:
        #     np.random.shuffle(index)
        # gene_list = gene_list.iloc[index, :]
        # dnv_table = dnv_table.iloc[index, :]
        # gene_scores = gene_scores.iloc[index, :]

        # self.x_dim = [gene_scores.shape[1], dnv_table]

        # increase labeled data points
        if train_flag and warmup_flag:
            x_labeled_index = np.where(~np.isnan(label.iloc[:, 0].values))[0]
            x_unlabeled_index = np.where(np.isnan(label.iloc[:, 0].values))[0]
            repeat_times = np.ceil(x_unlabeled_index.shape[0]/x_labeled_index.shape[0])
            x_labeled_index = np.repeat(x_labeled_index, repeats=repeat_times)
            x_index = np.concatenate((x_labeled_index, x_unlabeled_index))
            x_in = x_in.iloc[x_index, :]
            x_out = x_out.iloc[x_index, :]
            x_var = x_var.iloc[x_index, :]
            label = label.iloc[x_index, :]
        # shuffle x
        index = np.arange(x_in.shape[0])
        if train_flag:
            np.random.seed(0)
            np.random.shuffle(index)
        x_in = x_in.iloc[index, :]
        x_out = x_out.iloc[index, :]
        x_var = x_var.iloc[index, :]
        label = label.iloc[index, :]

        self.x_in = x_in
        self.x_out = x_out
        self.x_dim = x_dim
        self.x_var = x_var
        self.label = label
        self.train_flag = train_flag
        self.warmup_flag = warmup_flag

    def __iter__(self):
        """
        Iterates over all genes in the gene_list: self.gene_list in the data,
        yielding:
            x: node feature, in this case gene counts, with shape [num_nodes, 1]
            edge_index: single cell Graph connectivity in COO format with shape [2, num_edges]
            edge_attr: single cell Graph connectivity in COO format with shape [num_edges, 1]
            y: Graph classification
        """
        for i in range(self.x_in.shape[0]):
            yield (self.x_in.iloc[i, :].values.astype(np.float32),
                   self.x_out.iloc[i, :].values.astype(np.float32),
                   self.x_var.iloc[i, :].values.astype(np.float32),
                   self.label.iloc[i, :].values.astype(np.float32)
                   if self.train_flag and self.warmup_flag else np.array([np.nan]))

    def __len__(self):
        return len(self.x.shape[0])

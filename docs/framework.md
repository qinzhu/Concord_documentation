# The CONCORD Framework

## What Does CONCORD Do?

CONCORD addresses three foundamental challenges in single-cell data analysis — dimensionality reduction, denoising, and data integration—within a single unified framework. Operating in a fully unsupervised manner, it generates denoised, high-resolution, and batch-corrected embeddings that faithfully capture the structure of the underlying cell state landscape.

![State manifold](images/sub1.png){width="60%"}

In addition to its core functionality, CONCORD supports a range of downstream tasks, including cell type classification, doublet detection, cross-dataset projection, and annotation-guided representation learning. It also provides tools for simulating single-cell data and benchmarking dimensionality reduction and integration methods.

## What Makes CONCORD Powerful?

CONCORD's core innovation lies in its **mini-batch sampling framework**, which introduces a novel probabilistic strategy for sampling cells into mini-batches — the basic units of machine learning training. Combined with contrastive learning, this approach enables CONCORD to achieve state-of-the-art performance without relying on deep architectures, auxiliary losses, or supervision.

- **Neighborhood-Aware Sampling:**  
  Unlike conventional samplers, the neighborhood-aware sampler allows the model to explore local regions of the landscape while maintaining a global perspective. This enables the model to learn subtle distinctions among closely related cell states.

    ![Neighborhood-aware sampler](images/sub2.png){width="90%"}

- **Dataset-Aware Sampling:**  
  When applied to a single dataset, contrastive learning effectively captures biological variation in the latent space: 
  ![Single dataset contrastive learning](images/sub_singd.png){width="90%"}

    However, with uniform sampling across multiple datasets, both biological and dataset-specific variations are encoded, leading to latent spaces that separate by dataset as well as cell type:
    ![Multi dataset contrastive learning](images/sub_multid.png){width="90%"}
  
    To address this, we also introduce a dataset-aware sampler that restricts mini-batches to a single dataset, ensuring contrasts reflect only biological differences, as in the single-dataset setting:
    ![Multi dataset CONCORD](images/sub_intmulti.png){width="90%"}

    Dataset-specific biases are further diminished through random mini-batch shuffling: if such signals are encoded in one batch, they are disrupted and overwritten by subsequent mini-batches from other datasets. Consequently, only biologically meaningful signals, such as gene co-expression patterns, persist throughout training, resulting in a latent space that reflects biological variation with minimal batch effects. Importantly, this strategy of removing batch effect imposes no assumptions on the structure of the data beyond the existence of shared biological programs across datasets.

- **Joint probablistic Sampling:**  
  Both the neighborhood-aware and dataset-aware samplers follow a unified principle: probabilistically structuring mini-batches to balance global biological diversity with local and dataset-specific variation. We integrate both samplers into a joint sampling framework, where the likelihood of selecting a cell satisfies both sampling schemes:

    ![CONCORD sampler](images/sub4.png){width="90%"}

## The CONCORD Framework

The CONCORD framework is compatible with custom model architectures ([see API](../api/model/#concord.model.ConcordSampler)).  
In our study, we adopted a minimalist model: a single hidden-layer encoder with optional decoder and classifier heads.

![CONCORD model](images/sub5.png){width="70%"}

This minimalistic, highly extensible framework scales from small to large datasets and generalizes to modalities beyond scRNA-seq.

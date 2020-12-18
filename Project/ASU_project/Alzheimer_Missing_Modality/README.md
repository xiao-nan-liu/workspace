This document is to describe the Alzheimer's Disease Project.

Background

Use the MRI, FDG-PET, and Amyloid-PET image to detect if the patient is in earlystage of AD, which is, cognitive impairment.

Challenge

Not all patients have all these three types of images

Proposed Method

Use Incomplete Modality Transfer Learning approach to deal with the problem. More specifically, given the assumption that all patients have MRI image, we proposed the conditional distribution among different image features, and then use EM algorithm to handle missing modality problem to estimate parameters. See detailsin the Paper.

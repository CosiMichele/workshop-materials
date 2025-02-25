# Handling Images & Videos pt. 2

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/Geospatial_Workshops/main/images/coco_labels.png" width="600">
</p>
<br>

---

>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:10pm: Welcome and Introduction to the topic
> - 2:10pm-2:30pm: Overview of tools
> - 2:30pm-end: Testing and executing tools

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Understanding of Convolutional Neural Networks
> - Executing a small Jupyter Notebook using JS2

<br>

---
---

## Overview

<br>
<br>
<p align="center">
    <img src="https://camo.githubusercontent.com/31bd0e2e5b3a6e6abb8eaa6cb2f4284eab6ca727c649a24942a6381ccdc702f7/68747470733a2f2f696d616765732e73717561726573706163652d63646e2e636f6d2f636f6e74656e742f76312f3537663664353163396637343536366635356563663237312f64616564376631362d353237662d343135302d386264642d6362623230653236373435312f636865657461682d657a6769662e636f6d2d766964656f2d746f2d6769662d636f6e7665727465722e6769663f666f726d61743d31383077" width="400">
</p>
<br>

As discussed last week, Convolutional Neural Networks (CNNs) are a useful technique widely used to extrapolate information from images and videos. Today, we discuss some of its potential applications in science, such as object recognition and tracking. Additionally, we will also see how to set up your machine in order to execute a CNN model.

---

## Scientific Applications

Object tracking has been an active fields for decades, but with Open Science becoming more prominent and with hardware becoming more affordable, open source code is now much easier to obtain and apply. 

For example, in the field of biology there are at least 2 software available to the pubic: [DeepLabCut](https://github.com/DeepLabCut/DeepLabCut) and [SLEAP](https://github.com/talmolab/sleap).

### [DeepLabCut]((https://github.com/DeepLabCut/DeepLabCut))

<br>
<p align="center">
    <img src="https://camo.githubusercontent.com/73a9872de7208eba813b9668df10bc99ed1064294d5f8f6f56e025d33c9fecd8/68747470733a2f2f737461746963312e73717561726573706163652e636f6d2f7374617469632f3537663664353163396637343536366635356563663237312f742f3563336662656437346661353161636563643633646565622f313534373638313533343733362f4d6f7573654c6f636f6d6f74696f6e5f77617272656e2e6769663f666f726d61743d35303077" width="300">
</p>
<br>

DeepLabCut (DLC) is an open-source, deep-learning-based software designed for pose estimation. It enables researchers to track and analyze movements of animals in videos by detecting specific body parts without the need for physical markers.

It was inspired by the need for precise, automated tracking of animal behavior in neuroscience and behavioral research. The software is built upon [DeepLab](https://arxiv.org/abs/1606.00915) (Chen *et al.*, 2016), a deep-learning framework initially designed for image segmentation, utilizing ResNet (Residual Networks) as the backbone architecture. It employs transfer learning from pre-trained models on large image datasets, refining them with user-labeled frames for precise keypoint detection. DLC applies this framework to the problem of keypoint detection, allowing researchers to track body parts efficiently. [viso.ai has a quick rundown on what DeepLab is](https://viso.ai/deep-learning/deeplab/). 

DeepLabCut is widely used across various scientific fields, including:

- **Neuroscience** – Tracking animal movements to study motor functions and behavioral responses.
- **Biomechanics** – Analyzing motion and gait in animals and humans.
- **Ecology & Ethology** – Studying natural behaviors of animals in the wild or laboratory settings.
- **Medical & Rehabilitation Research** – Assessing motor impairments in preclinical models of neurological disorders.
- **Sports Science**– Analyzing human motion for performance improvement and injury prevention.
- **Human-Computer Interaction (HCI)** – Tracking human hand and facial movements for gesture recognition.

### [SLEAP](https://github.com/talmolab/sleap)

<br>
<p align="center">
    <img src="https://sleap.ai/docs/_static/sleap_movie.gif" width="400">
</p>
<br>

SLEAP (Social LEAP Estimates Animal Poses) is an open-source deep-learning-based tool for multi-animal pose estimation. It is designed as an alternative to DeepLabCut (DLC) but focuses on handling multiple interacting animals more effectively.

In order to help tracking multiple objects, SLEAP uses a similar setup as DLC but also includes multiple network architectures for different use cases:

- **Top-down models**: first detecting the entire animal and then estimating its pose.
- **Bottom-up models**: directly detecting individual keypoints across all animals.
- Fully convolutional networks for segmentation-based approaches.
- **ResNet** and **HRNet** architectures for feature extraction and pose estimation.

---

## The Swiss Army Knife (Resources for the Bigger Picture)

This workshop covers the tools required for the FOSS-AI goal of creating a small live object recognition demo. Here are the libraries required (and additional tools you may want to track):

### [YOLO](https://arxiv.org/abs/1506.02640)

*Read the paper*: [You Only Look Once](https://arxiv.org/abs/1506.02640), Redmond *et al.*, 2015

<br>
<p align="center">
    <img src="https://media.datacamp.com/legacy/image/upload/v1664382694/YOLO_Architecture_from_the_original_paper_ff4e5383c0.png" width="500">
</p>
<br>

YOLO (You Only Look Once) is a detection algorithm known for its speed and accuracy. Unlike traditional object detection methods that use sliding windows or region proposals, YOLO processes an image in a single forward pass, making it significantly faster. The goal was to create an object detection model that could work in real-time, unlike region-based models such as R-CNN or Faster R-CNN, which were slower.

Quoting [DataCamp](https://www.datacamp.com/blog/yolo-object-detection-explained): "*The authors frame the object detection problem as a regression rather than a classification task by spatially separating bounding boxes and associating probabilities to each detected image using a single convolutional neural network (CNN).*"

Additional Resource: [What is YOLO algorithm
, Medium](https://medium.com/@ishudey11032002/what-is-yolo-algorithm-ef5a3326510b).

### [Ultralytics](https://docs.ultralytics.com/quickstart/#install-ultralytics)

<br>
<p align="center">
    <img src="https://cdn.prod.website-files.com/646dd1f1a3703e451ba81ecc/64994690ba750f8de7df4478_Ultralytics_full_blue.svg" width="300">
</p>
<br>

Ultralytics is a company that develops and maintains the YOLO algorithms. We are going to be using Utralytics in order to train and run the YOLO models, as scripts and commands they have created do the heavy lifting.

One can install Ultralytics on their machine using

```
pip install ultralytics
```

or through conda

```
conda install -c conda-forge ultralytics
```

>[!important]
> :heavy_exclamation_mark: **Important!**
> For YOLO and Ultralytics to work as intended, one should have **[PyTorch](https://pytorch.org/)** installed.

### [Gradio](https://www.gradio.app/) 

<br>
<p align="center">
    <img src="https://pypi-camo.freetls.fastly.net/a95ef5913dc4cc84d2155ff690a0fa0d4c33d7e2/68747470733a2f2f7261772e67697468756275736572636f6e74656e742e636f6d2f67726164696f2d6170702f67726164696f2f6d61696e2f726561646d655f66696c65732f67726164696f2e737667" width="300">
</p>
<br>

Python library facilitating creating web apps to interface with an AI model.

Install with

```
pip install --upgrade gradio
```

### [OpenCV](https://pypi.org/project/opencv-python/)

<br>
<p align="center">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/d/d2/OpenCV_logo_black.svg/1200px-OpenCV_logo_black.svg.png" width="100">
</p>
<br>

OpenCV (Open Source Computer Vision Library) is a high-performance computer vision and image processing library. It provides tools for image manipulation, video analysis, object detection, and feature extraction.

In our case, it helps with creating visual output from interfacing with an AI model. Generates videos and images. Can also support webcam input for real-time AI.

```
pip install opencv-python
```

---

## Annotations

Annotations are essential for training machine learning models. Here are two resources that you can use for annotating a dataset.

### [Roboflow](https://app.roboflow.com/) 

<br>
<p align="center">
    <img src="https://d7umqicpi7263.cloudfront.net/img/product/8305253e-2066-4396-9e9a-f0f9b97e75b9.png" width="300">
</p>
<br>

Roboflow is a computer vision platform that helps users annotate, preprocess, and train object detection and classification models without requiring deep technical expertise. It's a website allowing image annotation, you can export the annotations for AI training for a number of formats (e.g., YOLO). Proprietary and requires a free account.

### [Label Studio](https://labelstud.io/) 

<br>
<p align="center">
    <img src="https://user-images.githubusercontent.com/12534576/192582340-4c9e4401-1fe6-4dbb-95bb-fdbba5493f61.png" width="400">
</p>
<br>

Open source image annotation, can also export image annotations in various formats.

---

**Additional Resources**:

- segment anything: https://segment-anything.com/
- [Follow this quick tutorial on running segment anything on your computer](https://www.hackster.io/lurst811/realtime-language-segment-anything-on-jetson-orin-ccf6e1)

---
---
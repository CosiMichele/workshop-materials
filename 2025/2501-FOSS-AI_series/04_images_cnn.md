# Handling Images & Videos pt. 1


📓 Notebook link: [CNN.ipynb](https://github.com/ua-datalab/MLWorkshops/blob/main/Convolutional_Neural_Networks/CNN.ipynb) ([raw](https://github.com/ua-datalab/MLWorkshops/blob/main/Convolutional_Neural_Networks/CNN.ipynb?raw=true))

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
> - 2:10pm-2:20pm: Convolutional Neural Network Architecture (+ *what's a Convolution?*)
> - 2:20pm-3:00pm: 

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

## Computer Vision Tasks

Generally, we want to automate the derivation of useful information from images. Some example tasks include:
- **Classification**: given an image, predict a class label
- **Object detection**: generate a bounding box around the object
- **Semantic segmentation**: assign every pixel in the image a class label
- **Instance segmentation**: differentiate between multiple instances of the same semantic class
- **Pose recognition**: for example, estimating the pose of a head, which can be used to determine what they are looking at
- **Activity recognition**: related to pose recognition, classify a pose or series of poses
- **Object tracking**: propose correspondence of detected objects across frames of a video
- **Image restoration**
- **Feature matching**: detection of features and correspondence between multiple views

---

## CNN Architecture

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/image.webp" width="350">
</p>
<br>

Traditional image analysis approaches such as looking at color values of individual pixels or grouping adjacent pixels based on similar values (i.e., objects), will both likely to fail at the task of identifying objects in high-resolution images.

A convolutional neural network (CNNs) can be trained to identify objects in an image more analogously to how humans can identify objects. We generally look at multiple features of an object to decide what the object is. For example, if presented with an image of a tree, we may look at features such as crown shape, leaf architecture, and color to help us identify the object as a tree.

Please watch these VERY useful videos to understand how Convolutional Neural Networks operate:

- [*What are Convolutional Neural Networks (CNNs)?*](https://www.youtube.com/watch?v=QzY57FaENXg) from IBM Technology
- [*But what is a neural network?*](https://www.youtube.com/watch?v=aircAruvnKk&list=PLZHQObOWTQDNU6R1_67000Dx_ZCJB-3pi) from 3Blue1Brown (*highly recommended!!*)
- [*How Convolutional Neural Networks Work*](https://www.youtube.com/watch?v=FmpDIaiMIeA) from Brandon Rohrer
- [*Convolutional Neural Networks Explained (CNN Visualized)*](https://www.youtube.com/watch?v=pj9-rr1wDhM) from Futurology

In this workshop, we will use RGB images. Each image has dimensionality H x W x 3, where H is the height of the image, W is the width of the image, and every pixel has three color channels (Red, Green, Blue).

To motivate the need for a neural network architecture for images, consider the case of classifying digits with MNIST using a fully connected neural network. In this case each pixel takes on only a single greyscale value, and thus the image has H x W values. In the fully connected neural network, we unravel the image into a one dimensional vector which is HW long, either by picking values row-wise or column-wise. 

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/vector.JPG" width="350">
</p>
<br>


Despite both images having the same structure, the location of the white value is shifted in the vector, meaning it interacts with a completely different set of weights and biases. This is a toy example of a larger problem. For a classification example, we would like to learn to identify an image that features a cat, whether the cat is in the upper left or the lower right of the image. This property is known as [translation invariance](https://stats.stackexchange.com/questions/208936/what-is-translation-invariance-in-computer-vision-and-convolutional-neural-netwo). Convolution is the linear operator that enables this in CNNs.

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/architecture_overview.webp" width="600">
</p>
<br>


Convolutional neural network architectures consist of convolutional and pooling layers in addition to fully connected layers. Note that we use a common activation function, ReLU. Based on the output of the fully connected layers, we are performing a classification task between some number of classes. Note also that the network reduces the dimensionality of the input using the convolution and pooling layers to some compact representation, and then feeds this into the fully connected layers which will generate the predicted class.

---

### Convolution

In the convolution operation, we pass kernels (also known as filters) over an image. We take the inner product between the kernel and a patch of the image. The output of convolution is high where the underlying image resembles the filter, and low where it does not. Once this response is calculated, it is fed through an activation function like ReLU, similar to fully connected neural networks.

<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/kernel_stationary.png" width="250">
</p>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/kernel_sweep.gif" width="250">
</p>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/kernel_sweep_3d.webp" width="250">
</p>
<br>

Before neural networks there was significant effort expended to design kernels to detect small features, such as edges, corners, etc and use these for computer vision tasks. Neural networks allow us to set up a convolutional architecture that both learns the kernels that are useful for the given task as well as the mapping from the feature space to the output - i.e. the network learns the parameters of each kernel.

Here's an example of the learned filter bank of ImageNet (Krizhevsky et al.). Note that some of the filters are what we would expect: the network has learned to look for lines at various angles, as well as for dots. Not all of the filters are easily interpretable.

![](https://github.com/ua-datalab/MLWorkshops/blob/main/images/CNN/ImageNet_kernels.jpeg)

In addition to specifying the number of filters and the shape of the filters, we also need to specify the stride and padding. The stride specifies how many pixels the kernel shifts by as it slides around the input image. Using zero-padding, the image will add artificial black pixels in a border around the input.

![](https://github.com/ua-datalab/MLWorkshops/blob/main/images/CNN/stride_padding.png)

---

### Pooling

These feature maps still have very high dimension, so we pool the filter activations through a pooling layer. This can be thought of as summarizing the input. In max pool, we take the maximum of all activations, and in average pool we take the average of all activations. Typically max pool is used. There are no learnable parameters for this layer.

<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/pool.jpeg" width="250">
</p>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/pooling_types.webp" width="250">
</p>
<br>

---

### Parameter Sharing

In the AlexNet architecture, which represented a huge leap forward in image classification, input images are 227x227x3. The first convolutional layer output is 55x55x96 (290,400), so 96 filters must be learned. Each filter is 11x11x3, so we learn 11x11x3 weights and 1 bias for a total of 364 parameters per filter. Thus we need to learn 34,944 parameters in total for this convolutional layer.

If we used a fully connected neural network instead, we would go from 227x227x3 (154,587) -> 55x55x96 (290,400), which requires 154587*290400 weights and 290400 biases. That's approximately 44 billion parameters for a single layer!

Also note that by learning filters that are shared across all spatial locations in the image, we have achieved our goal of translation invariance.

---

### Activations

<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/convnet_car_activations.jpeg" width="500">
</p>
<br>

Intuition: the first layers of the network learn low level, simple features. As the input progresses deeper into the architecture of the network, these are combined to create more complex features and eventually a high level semantic description of the image is extracted. This high level semantic description is then used in the fully connected network for the learning task.

## Case study: VGG-16 Architecture

On the ImageNet visual database, can classify images into 1000 classes with 92.7% accuracy. There are better models now, which use more advanced architectures, but this is quite simple and effective.

<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/MLWorkshops/refs/heads/main/images/CNN/vgg-16-cnn-architecture.webp" width="500">
</p>
<br>

## Sources for images and additional resources

- DataLab Machine Learning Workshop series (2024): [Convoultional Neural Networks](https://github.com/ua-datalab/MLWorkshops/wiki/Convolutional-Neural-Networks)
- DataLab Geospatial Workshop series (2024): Image Object Detection ‐ [Deep Forest](https://github.com/ua-datalab/Geospatial_Workshops/wiki/Image-Object-Detection-%E2%80%90-Deep-Forest) and [Detecto](https://github.com/ua-datalab/Geospatial_Workshops/wiki/Image-Object-Detection-%E2%80%90-Detecto)
- [Stanford's CS231n, CNNs for Visual Recognition](https://cs231n.github.io/) An excellent resource for neural networks as well as CNNs
- [ConvNetJS CIFAR-10](https://cs.stanford.edu/people/karpathy/convnetjs/demo/cifar10.html) Train and test a model on CIFAR10 in your browser. Excellent demonstration of the activations at every step in the architecture.
- [Comprehensive Guide to Convolutional Neural Networks - the ELI5 way](https://towardsdatascience.com/a-comprehensive-guide-to-convolutional-neural-networks-the-eli5-way-3bd2b1164a53)
- [A guide to convolution arithmetic for deep learning](https://arxiv.org/abs/1603.07285) A thorough overview of convolution, stride, and padding
- [Introduction to CNN](https://github.com/nextgensh/ua-teaching/blob/main/cnn/cnn-intro.ipynb). Shravan Aras.
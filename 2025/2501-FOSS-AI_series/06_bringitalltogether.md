# Bringing It All Together

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/ua-datalab/Geospatial_Workshops/main/images/coco_labels.png" width="600">
</p>
<br>

---

>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:10pm: 
> - 2:10pm-2:30pm: 
> - 2:30pm-end: 

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Refresh required Open Science libraries for image detection
> - Together train and deploy a model 

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

Today is the day we bring it all together! Over the course of this and the next workshop (if needed) we are going to take all of the Open Science pieces we have discussed through the [Functional Open Science Skills for AI/ML Applications](https://github.com/ua-datalab/FOSS_AI-ML/wiki) and create and end-to-end ML pipeline that will recognize objects within a photo or a video.

As a summary, here are the contents that we have covered and we are going to be using today:

- [Virtual Environments](https://docs.python.org/3/library/venv.html) (through Python)
- [Where to find compute](https://docs.jetstream-cloud.org/) (through JetStream2)
- Open Source ML libraries ([YOLO](https://docs.ultralytics.com/models/yolo11/)/[Ultralytics](https://docs.ultralytics.com/))
- Tools to build annotations for training, testing and deploying an ML object detector (through [Roboflow](https://roboflow.com/) (requires account)) 

Here's the general overview of what we are going to be doing today:

```mermaid
graph LR;
    'Accessing JetStream2 VMs'-->'Creating a virtual evironment';
    'Creating a virtual evironment'-->Installing required packages';
    'Installing required packages'-->'Training YOLO';
    'Annotate/Access Annotations (Roboflow)'-->'Training YOLO';
    'Traning YOLO'-->'Deploy trained weights on image/video';
```

---
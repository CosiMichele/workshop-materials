Sources:

- https://github.com/ua-datalab/MLWorkshops/wiki/Convolutional-Neural-Networks
- https://github.com/ua-datalab/Geospatial_Workshops/wiki/Image-Object-Detection-%E2%80%90-Detecto
- https://www.v7labs.com/blog/yolo-object-detection
- https://www.datacamp.com/blog/yolo-object-detection-explained

Tasks:
 
- Create image on JS2 with cv2, ultralytics
- Gather images of beans

Extra resources:

- Python Libraries:
  - [ultralytics](https://docs.ultralytics.com/quickstart/#install-ultralytics) - used for training and running yolo models, does most of the heavy lifting and can be used in scripts or commands
  - [gradio](https://www.gradio.app/) - Python library facilitating creating web apps to interface with an AI model
  - [opencv (denoted cv2 when importing)](https://pypi.org/project/opencv-python/) - Helps with creating visual output from interfacing with an AI model. Generates videos and images. Can also support webcam input for real-time AI.
- Docker Containers:
  - [jetson-inference](https://github.com/dusty-nv/jetson-inference) - Docker with python scripts, JetPack (python ai libraries that work on jetson (i.e. pytorch). This container is kind of tricky since it's jetson exclusive and its only real value comes from training classification models with pre-written python scripts
- Image annotation:
  - [Roboflow](https://app.roboflow.com/) - Website allowing image annotation, you can export the annotations for AI training for a number of formats, I only used yolo though. Proprietary and requires a free account
  - [Label Studio](https://labelstud.io/) - Open source image annotation, can also export image annotations in various formats. Not as familiar with it myself, however.

  Plan:

  - Day 1: discuss the general approaches. Images, CNN. Example of easy to find resources. Run example on JS2, using Jeff code.
  - Day 2: talk about videos, run example on JS2, use Leonardo's code.
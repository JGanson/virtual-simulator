# Project Introduction
The aim of this project is to create a virtual environment and leverage deep learning techniques to train a model capable of solving tasks within that environment. The project involves the integration of a virtual environment, possibly using OpenGL, where the simulation of a robot-like entity can take place. This environment will serve as the testing ground for training a deep learning model.

The focus of the project is to utilize reinforcement learning, a subfield of deep learning, to train the model. The reinforcement learning paradigm enables the model to learn from its interactions with the virtual environment, receiving rewards or penalties based on its performance in solving the designated task.

By training the deep learning model within the virtual environment, the project aims to enhance the model's ability to tackle specific tasks, such as walking to a designated position. The training process will involve iterations of learning, where the model progressively improves its performance through trial and error, optimizing its actions based on the provided rewards or penalties.

The successful implementation of this project will demonstrate the effectiveness of integrating a virtual environment and deep learning techniques. It has the potential to contribute to the development of intelligent systems capable of autonomously solving tasks within virtual environments, opening avenues for applications in fields such as robotics, artificial intelligence, and simulation.

<br></br>
# Simulation Environment
We will use C++ OpenGL for this part.

## Main Goals: 
<li> Be able to load and editing mesh or object model
<li> Have collision detection (Not necessary)

## Further Implementation
<li> Ray Tracing
<li> Force System Implementation


<br></br>
# Deep Learning
We will use Pytorch to implement our deep learning model. 
## Model
<li> Maybe a CNN or CV-Transformer to do vision detection
<li> Other Layers to let model perform different action based on the CV input

## Learning
<li> Finding a suitable reward function
<li> Finding a sutable Reinforcement learning algorithm 
<li> Understanding before implemented!

<br></br>
# Task For Summer
Our primary task for the summer is to load a car into the virtual environment along with a flat ground plane. We will designate a specific region on the ground with a distinct color and instruct the car to navigate to that region. Notably, the car should only rely on camera information for movement, without any additional data such as environmental position coordinates.

Additionally, we can explore more complex tasks, such as guiding the robot to move along a colored line on the ground.
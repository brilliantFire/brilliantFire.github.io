---
layout: post
comments: true
tags: TensorFlow GPU Anaconda Windows10
title: Installing TensorFlow-GPU on a Windows 10 Machine
---
A few months back I discovered that my three-year-old Dell desktop has a small but nonetheless CUDA-enabled NVIDIA GPU. While the [GeForce GT 730](https://www.geforce.com/hardware/desktop-gpus/geforce-gt-730) I found - with it's itty-bitty set of 384 cores and tiny 2GB stash of VRAM -  is not super-great for training huge, gnarly deep neural networks, I figured any GPU is better than no GPU. So, as an aspiring data scientist and big, big fan of neural networks (both real and artificial), I got a little excited, wondering... 

> #### "How do I get TensorFlow running on this bad Larry *RIGHT EFFING NOW!?!?*"

Turns out, this is not at all a clear-cut endeavor if your operating system is Windows 10. I did get it working but "RIGHT EFFING NOW" ended up being several hours of Googling, downloading, installing, swearing, uninstalling, swearing some more, and then, finally, finding some luck. The process that eventually worked for me (details below) was pieced together from various posts ([here's one](http://www.netinstructions.com/how-to-install-and-run-gpu-enabled-tensorflow-on-windows/)) and [the official instructions for installing the CUDA Toolkit on Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/). Your mileage may vary (read: I'm sure this won't be everyone's solution).

## Some Preliminaries

#### 1. Check out [this list](https://developer.nvidia.com/cuda-gpus) and verify that you have a CUDA-enabled GPU with a "Compute Compatability" of *3.0 or higher*.
You can install and use TensorFlow-GPU on a machine with a GPU of lower compute compatability just fine except that it won't actually run things on the GPU. Just the CPU. 

#### 2. The installation I'm about to describe here requires [Anaconda](https://www.anaconda.com/download/), a package and environment manager for Python. 
Click the link, then download and install the Python 3.6 version, if you don't already have it. I love Anaconda and could say a lot about it but I'll save that for another post <3

#### 3. (Optional) I recommend downloading [Revo Uninstaller](https://www.revouninstaller.com/revo_uninstaller_free_download.html?gclid=CjwKCAjwnLjVBRAdEiwAKSGPI481vBfbg8ZZL9fXGXW5v4c9Zbnk8Y-KJ0mUoeBcdYtLvde_4AxyZhoCMnEQAvD_BwE) and using it to uninstall CUDA if you have to go back and try again.
If the CUDA installation doesn't work for some reason, I found it to be easiest to uninstall it and try again. Using something like Revo Uninstaller to clear everything out (directories, etc.) helps.  

## 7 Steps to Installing TensorFlow-GPU

#### 1. Download and install [Microsoft Visual Studio Community 2015 with Update 3](https://my.visualstudio.com/Downloads?q=Visual%20Studio%202015%20with%20Update%203).
The link above will take you to a page that requires a free developer account with Microsoft. Follow the recommended installation instructions. This could take a while.

#### 2. Download and install [CUDA 9.0](https://developer.nvidia.com/cuda-90-download-archive?target_os=Windows&target_arch=x86_64&target_version=10&target_type=exenetwork) (_NOT_ 9.1).
As of this post (April 2018), 9.1 is the most recent version of the CUDA Toolkit (but it looks like 9.2 is coming soon). TensorFlow-GPU (I have version 1.1.0) right now only works with CUDA 9.0, so be sure to get that version. Follow the recommended installation.

#### 3. Open Python (command line/IDE/JupyterLab/Notebook - doesn't matter how) and try to import tensorflow.
Open Python and enter `import tensorflow`. It's going to fail but, like most failures, it should teach us something. In this case,  
- Whether our CUDA installation worked; and, if it did,
- Which version of the cuDNN library we need to download and put into the CUDA installation.  

##### *If the errors say something about missing "cuDNN64_7.dll" files or something like that...*
CONGRATULATIONS!! CUDA appears to be successfully installed! What's more, we can tell which version of the cuDNN library we need by looking at the number after the underscore in the name of the file it's asking for. This should be 7, which is the cuDNN version noted in [TensorFlow's instructions for installing on Windows](https://www.tensorflow.org/install/install_windows). Keep this in mind and proceed to Step 4.  

##### *If the errors mention a missing "cudart64_90.dll" file or anything other than cuDNN...*
Well, we can either:
- Proceed to Step 4 and try manually setting the `Path` environment variable (YIKES), or
- Uninstall and reinstall CUDA until it works (my recommendation).

#### 4. Double-check the `Path` environment variable.
If everything went swimmingly for you in Step 3, these should be all set up for you. If not, you can set them up manually and it works (I tried it). Either way, it's a good idea to check them.

To find the environment variables, we want to open the System Properties menu (type "control sysdm.cpl" without the quotes into a shell, or search for "System" in Windows Settings and click the link for "Edit the system environment variables"). Click the Environment Variables button on the Advanced tab. We're looking for the system variable named `Path` (green box in the in image below). Click on it and then hit Edit (orange box).

![Environment Variables](https://i.imgur.com/fc7LHby.png)

You should see a few instances of the CUDA directory in there (green box below). If you don't, you can add them to `Path` by clicking "New".

![Path Settings](https://i.imgur.com/vUYB6ZG.png)

#### 5. Download the [cuDNN 7.0.5 (Dec 5, 2017 release) library](https://developer.nvidia.com/rdp/cudnn-archive).
You'll need to sign up for a free developer account with NVIDIA for access. The above link goes to the archive site. Scroll down and click the "Download cuDNN v7.0.5 (Dec 5, 2017), for CUDA 9.0" link. Download the Windows 10 version. It's a .zip file. Save it someplace memorable. Once the download is finished, unzip in the directory you saved it to.

#### 6. Add the cuDNN files to the CUDA installation.
Once you've unzipped the cuDNN download, you'll see that you have a "CUDA" folder with three more folders inside called 'bin', 'include', and 'lib'. We need to find the corresponding 'bin', 'include', and 'lib' folders in the CUDA installation which will be located somewhere like
`
C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v.9.0
`
Take the cuDNN file from the 'bin' folder (the one you got by unzipping) and drag/drop it into the 'bin' folder in the CUDA installation. Repeat for the cuDNN files in 'include' and 'lib'. **Do not drag and drop the entire folder.** Only the files inside.

#### 7. Create a conda environment with Python 3.5 and install the TensorFlow-GPU package.
From here, I mostly followed the instructions in the **Installing with Anaconda** section of the [Windows installation instructions for TensorFlow](https://www.tensorflow.org/install/install_windows). The only thing I changed was to add the [Anaconda suite of packages](https://docs.anaconda.com/anaconda/packages/py3.5_win-64) to the environment (that way things like pandas and sklearn are already in there).  

Open the Anaconda prompt and use the following to create the environment, named "tensorflow":
`
C:> conda create -n tensorflow pip python=3.5 anaconda
`
(The "anaconda" on the end adds those packages I was talking about.)

Next, activate the environment by entering:
`
C:> activate tensorflow
`
The prompt should change. Now, use the following to install the tensorflow-gpu package:
`
(tensorflow)C:> pip install --ignore-installed --upgrade tensorflow-gpu 
`

##### And you should be good to go!
__One last important tip:__ If you're going on to install the [Keras](https://keras.io/) package, be sure to install the *GPU version (keras-gpu)*. The other version comes with non-GPU TensorFlow as a dependency and that will overwrite the GPU version you just installed.

### Test the installation.
You can use this code to see whether TensorFlow is installed and that it's using your GPU.
`
from tensorflow.python.client import device_lib
device_lib.list_local_devices()
`
If it is, you should see your GPU in the list (like the GeForce GT 730 below).
![SUCCESS!](https://i.imgur.com/8AImXBR.png)

WOOHOO! 

![HIGH FIVE!](https://media.giphy.com/media/BwOU6uH7afefu/giphy.gif)

Happy network building!

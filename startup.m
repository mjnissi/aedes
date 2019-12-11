%% This is a script m-file that will be run each time Matlab is started,
%% provided that this file is in the Matlab search path (see PATH
%% command).

% If there appears any graphical anomalities when running AEDES, try to
% some of the following OpenGL workarounds. If you experience some of the
% symptoms descriped below uncomment the line under it and restart Matlab.

%% -----------------------------------------------------------------

% OpenGLClippedImageBug Symptom: Images (as well as colorbar displays) do
% not display when the Renderer property set to opengl. In AEDES this is
% seen as missing images, i.e. there are no images and no error messages
% when data is opened into AEDES. This seams to be very common issue with
% some integrated graphics found mostly in Notebook computers (especially
% SiS graphics adapters).

%opengl('OpenGLClippedImageBug',1)

%% -----------------------------------------------------------------

% OpenGLEraseModeBug
% Symptom: Erasemodes other than "normal" do not function properly. In
% AEDES this is seen as flashing of the text and crossbars over the
% image. This seems to be common with NVidia (GeForce) graphics
% adapters.

%opengl('OpenGLEraseModeBug',1)

%% -----------------------------------------------------------------

% OpenGLBitmapZbufferBug:
% Symptom: text with background color (including data tips) and text
% displayed on image, patch, or surface objects is not visible when using
% OpenGL renderer.

%opengl('OpenGLBitmapZbufferBug',1)

%% -----------------------------------------------------------------

% OpenGLWobbleTesselatorBug:
% Symptom: Rendering complex patch object causes segmentation violation
% and returns a tesselator error message in the stack trace. (Probably
% does not occur with AEDES, because it doesn't use complex patch objects)

%opengl('OpenGLWobbleTesselatorBug',1)

%% -----------------------------------------------------------------

% OpenGLLineSmoothingBug
% Symptom: Lines with a LineWidth greater than 3 look bad.

%opengl('OpenGLLineSmoothingBug',1)

%% -----------------------------------------------------------------

% OpenGLDockingBug
% Symptom: MATLAB crashes when you dock a figure that has its Renderer
% property set to opengl.

%opengl('OpenGLDockingBug',1)




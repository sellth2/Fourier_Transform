function varargout = four_tr(varargin)
% FOUR_TR MATLAB code for four_tr.fig
%      FOUR_TR, by itself, creates a new FOUR_TR or raises the existing
%      singleton*.
%
%      H = FOUR_TR returns the handle to a new FOUR_TR or the handle to
%      the existing singleton*.
%
%      FOUR_TR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOUR_TR.M with the given input arguments.
%
%      FOUR_TR('Property','Value',...) creates a new FOUR_TR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before four_tr_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to four_tr_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help four_tr

% Last Modified by GUIDE v2.5 22-May-2022 13:48:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @four_tr_OpeningFcn, ...
                   'gui_OutputFcn',  @four_tr_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before four_tr is made visible.
function four_tr_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to four_tr (see VARARGIN)

% Choose default command line output for four_tr
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes four_tr wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = four_tr_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value=get(handles.slider1,'value');
string_value=num2str(value);
set(handles.edit1,'String',string_value);
handles.o=value;
guidata(hObject, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in ideal_filter.
function ideal_filter_Callback(hObject, eventdata, handles)
% hObject    handle to ideal_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of ideal_filter
m=220;   % Size is set
n=220;
d0=handles.o;
a=im2double(rgb2gray(imread('lena.png')));
A = fft2(a);                                        %fourier transform of image
A1 = fftshift(A);                                   %shifting origin
Aabs = abs(A1);                                     %Magnitude of A1 (Frequency domain representation of image)

if get(handles.low_pass,'Value')
       
    for i=1:m
        for j=1:n
            d(i,j) = sqrt((i-m/2)^2+(j-n/2)^2);
            if d(i,j) < d0
                 Bl(i,j) = A1(i,j);
                 lowpass(i,j) = 1;
            else
                    Bl(i,j) = 0;
                    lowpass(i,j) = 0;

            end
       end
    end
    B1l = fftshift(Bl);                               %Reshifting the origin of filtered image
    bl = ifft2(B1l);                                  %Taking inverse fourier transform
    babs = abs(bl);                                   %Taking magnitude. (low pass Filtered output image)

    axes(handles.axes2);
    imshow(lowpass);
    handles.J=lowpass;
    guidata(hObject, handles);

    axes(handles.axes3);
    surf(lowpass);
    handles.J=lowpass;
    guidata(hObject, handles)

    axes(handles.axes4);
    imshow(babs);
    handles.B=babs;
    guidata(hObject, handles);

else
    for i=1:m
        for j=1:n
            dh(i,j) = sqrt((i-m/2)^2+(j-n/2)^2);
            if dh(i,j) > d0
                Bh(i,j) = A1(i,j);
                highpass(i,j) = 1;
            else
                Bh(i,j) = 0;
                highpass(i,j) = 0;

            end
        end
    end
    B1h = fftshift(Bh);                               %Reshifting the origin of filtered image
    bh = ifft2(B1h);                                  %Taking inverse fourier transform
    babsh = abs(bh);                                   %Taking magnitude. (High pass Filtered output image)

    axes(handles.axes2);
    imshow(highpass);
    handles.J=highpass;
    guidata(hObject, handles);

    axes(handles.axes3);
    surf(highpass);
    handles.J=highpass;
    guidata(hObject, handles)

    axes(handles.axes4);
    imshow(babsh);
    handles.B=babsh;
    guidata(hObject, handles);
end

% --- Executes on button press in gaussian_filter.
function gaussian_filter_Callback(hObject, eventdata, handles)
% hObject    handle to gaussian_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m=220;   % Size is set
n=220;
d0=handles.o;
a=im2double(rgb2gray(imread('lena.png')));
A = fft2(a);                                        %fourier transform of image
A1 = fftshift(A);                                   %shifting origin
Aabs = abs(A1);                                      %Magnitude of A1 (Frequency domain representation of image)
if get(handles.low_pass,'Value')
    for i=1:m
        for j=1:n
            dgauss(i,j) = (exp(-(i-m/2)^2/(2*d0^2))*exp(-(j-n/2)^2/(2*d0^2)));
        end
    end

    Bgaussl = dgauss.*A1;
    Bgausslmag = dgauss.*Aabs;
    Bgaussl1 = fftshift(Bgaussl);
    bgaussl = ifft2(Bgaussl1);

    axes(handles.axes2);
    imshow(dgauss);
    handles.J=dgauss;
    guidata(hObject, handles);

    axes(handles.axes3);
    surf(dgauss);
    handles.J=dgauss;
    guidata(hObject, handles)

    axes(handles.axes4);
    imshow(bgaussl);
    handles.J=bgaussl;
    guidata(hObject, handles);

else
      % Gaussian High pass filtering

    for i=1:m
        for j=1:n
            dgaussh(i,j) = -(exp(-(i-m/2)^2/(2*d0^2))*exp(-(j-n/2)^2/(2*d0^2)))+1;
        end
    end

    Bgaussh = dgaussh.*A1;
    Bgaussh1 = fftshift(Bgaussh);
    bgaussh = ifft2(Bgaussh1);

    axes(handles.axes2);
    imshow(dgaussh);
    handles.J=dgaussh;
    guidata(hObject, handles);

    axes(handles.axes3);
    surf(dgaussh);
    handles.J=dgaussh;
    guidata(hObject, handles)

    axes(handles.axes4);
    imshow(bgaussh);
    handles.J=bgaussh;
    guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of gaussian_filter
end
% --- Executes on button press in butter_filter.
function butter_filter_Callback(hObject, eventdata, handles)
% hObject    handle to butter_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m=220;   % Size is set
n=220;
d0=handles.o;
a=im2double(rgb2gray(imread('lena.png')));
N=1;                        % Set the filter order
% Determine the centre of the image
p= round(m/2);
q= round(n/2);
if get(handles.low_pass,'Value')
    % Define the filter kernel
    H= zeros(m,n);
    for i = 1:m
        for j = 1:n
            d = (i-p).^2+(j-q).^2;
            H(i,j) = 1/(1+((d/d0/d0).^(2*N)));
        end
    end
    %Input Image in frequency domain
    A_f = fftshift(fft2(a));
    
    %Apply the filter
    B = A_f.*H;
    C = abs(ifft2(B));

    axes(handles.axes2);
    imshow(H);
    handles.J=H;
    guidata(hObject, handles);

    axes(handles.axes3);
    surf(H);
    handles.J=H;
    guidata(hObject, handles);

    axes(handles.axes4);
    imshow(C);
    handles.J=C;
    guidata(hObject, handles);


else
    %Define the kernel
    H= zeros(m,n);
    for i = 1:m
        for j = 1:n
            d = (i-p).^2+(j-q).^2;
            if d~=0
               H(i,j) = 1/(1+((d0*d0/d).^(2*N)));
            end
        end
    end
    %Input image in frequency domain
    A_f = fftshift(fft2(a));
    %Apply the filter
    B = A_f.*H;
    C = abs(ifft2(B));
    %Display
    axes(handles.axes2);
    imshow(H);
    handles.J=H;
    guidata(hObject, handles);

    axes(handles.axes3);
    surf(H);
    handles.J=H;
    guidata(hObject, handles);

    axes(handles.axes4);
    imshow(C);
    handles.J=C;
    guidata(hObject, handles);

end
% Hint: get(hObject,'Value') returns toggle state of butter_filter

% --- Executes on button press in low_pass.
function low_pass_Callback(hObject, eventdata, handles)
% hObject    handle to low_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get(handles.low_pass,'Value')
        set(handles.high_pass,'Value',0)
        if get(handles.ideal_filter,'Value')
            m=220;   % Size is set
            n=220;
            d0=handles.o;
            a=im2double(rgb2gray(imread('lena.png')));
            A = fft2(a);                                        %fourier transform of image
            A1 = fftshift(A);                                   %shifting origin
            Aabs = abs(A1);                                     %Magnitude of A1 (Frequency domain representation of image)
            for i=1:m
                for j=1:n
                    d(i,j) = sqrt((i-m/2)^2+(j-n/2)^2);
                    if d(i,j) < d0
                        Bl(i,j) = A1(i,j);
                        lowpass(i,j) = 1;
                    else
                        Bl(i,j) = 0;
                        lowpass(i,j) = 0;

                    end
                end
            end
            B1l = fftshift(Bl);                               %Reshifting the origin of filtered image
            bl = ifft2(B1l);                                  %Taking inverse fourier transform
            babs = abs(bl);                                   %Taking magnitude. (low pass Filtered output image)

            axes(handles.axes2);
            imshow(lowpass);
            handles.J=lowpass;
            guidata(hObject, handles);

            axes(handles.axes3);
            surf(lowpass);
            handles.J=lowpass;
            guidata(hObject, handles)

            axes(handles.axes4);
            imshow(babs);
            handles.B=babs;
            guidata(hObject, handles);

        elseif get(handles.gaussian_filter,'Value')
              m=220;   % Size is set
              n=220;
              d0=handles.o;
              a=im2double(rgb2gray(imread('lena.png')));
              A = fft2(a);                                        %fourier transform of image
              A1 = fftshift(A);                                   %shifting origin
              Aabs = abs(A1);                                      %Magnitude of A1 (Frequency domain representation of image) 
              for i=1:m
                  for j=1:n
                      dgauss(i,j) = (exp(-(i-m/2)^2/(2*d0^2))*exp(-(j-n/2)^2/(2*d0^2)));
                  end
              end

              Bgaussl = dgauss.*A1;
              Bgausslmag = dgauss.*Aabs;
              Bgaussl1 = fftshift(Bgaussl);
              bgaussl = ifft2(Bgaussl1);

              axes(handles.axes2);
              imshow(dgauss);
              handles.J=dgauss;
              guidata(hObject, handles);

              axes(handles.axes3);
              surf(dgauss);
              handles.J=dgauss;
              guidata(hObject, handles)

              axes(handles.axes4);
              imshow(bgaussl);
              handles.J=bgaussl;
              guidata(hObject, handles);
        else 
            m=220;   % Size is set
            n=220;
            d0=handles.o;
            a=im2double(rgb2gray(imread('lena.png')));
            N=1;                        % Set the filter order
            % Determine the centre of the image
            p= round(m/2);
            q= round(n/2);
            % Define the filter kernel
            H= zeros(m,n);
            for i = 1:m
                for j = 1:n
                    d = (i-p).^2+(j-q).^2;
                    H(i,j) = 1/(1+((d/d0/d0).^(2*N)));
                end
            end
            %Input Image in frequency domain
            A_f = fftshift(fft2(a));
    
            %Apply the filter
            B = A_f.*H;
            C = abs(ifft2(B));

            axes(handles.axes2);
            imshow(H);
            handles.J=H;
            guidata(hObject, handles);

            axes(handles.axes3);
            surf(H);
            handles.J=H;
            guidata(hObject, handles);

            axes(handles.axes4);
            imshow(C);
            handles.J=C;
            guidata(hObject, handles);
        end
 
  end
% Hint: get(hObject,'Value') returns toggle state of low_pass

% --- Executes on button press in show_image

% --- Executes during object creation, after setting all properties.
function low_pass_CreateFcn(hObject, eventdata, handles)
set(hObject,'Value',1); % if you want it to be seelcted initialy
% hObject    handle to low_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in high_pass.
function high_pass_Callback(hObject, eventdata, handles)
% hObject    handle to high_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.high_pass,'Value')
        set(handles.low_pass,'Value',0)
        if get(handles.ideal_filter,'Value')
            m=220;   % Size is set
            n=220;
            d0=handles.o;
            a=im2double(rgb2gray(imread('lena.png')));
            A = fft2(a);                                        %fourier transform of image
            A1 = fftshift(A);                                   %shifting origin
            Aabs = abs(A1);                                     %Magnitude of A1 (Frequency domain representation of image) 
            for i=1:m
                for j=1:n
                    dh(i,j) = sqrt((i-m/2)^2+(j-n/2)^2);
                    if dh(i,j) > d0
                        Bh(i,j) = A1(i,j);
                        highpass(i,j) = 1;
                    else
                        Bh(i,j) = 0;
                        highpass(i,j) = 0;

                    end
                end
            end
            B1h = fftshift(Bh);                               %Reshifting the origin of filtered image
            bh = ifft2(B1h);                                  %Taking inverse fourier transform
            babsh = abs(bh);                                   %Taking magnitude. (High pass Filtered output image)

            axes(handles.axes2);
            imshow(highpass);
            handles.J=highpass;
            guidata(hObject, handles);

            axes(handles.axes3);
            surf(highpass);
            handles.J=highpass;
            guidata(hObject, handles)

            axes(handles.axes4);
            imshow(babsh);
            handles.B=babsh;
            guidata(hObject, handles);
        elseif get(handles.gaussian_filter,'Value')
              m=220;   % Size is set
              n=220;
              d0=handles.o;
              a=im2double(rgb2gray(imread('lena.png')));
              A = fft2(a);                                        %fourier transform of image
              A1 = fftshift(A);                                   %shifting origin
              Aabs = abs(A1);                                      %Magnitude of A1 (Frequency domain representation of image) 
              % Gaussian High pass filtering

                for i=1:m
                    for j=1:n
                        dgaussh(i,j) = -(exp(-(i-m/2)^2/(2*d0^2))*exp(-(j-n/2)^2/(2*d0^2)))+1;
                    end
                end

                Bgaussh = dgaussh.*A1;
                Bgaussh1 = fftshift(Bgaussh);
                bgaussh = ifft2(Bgaussh1);

                axes(handles.axes2);
                imshow(dgaussh);
                handles.J=dgaussh;
                guidata(hObject, handles);

                axes(handles.axes3);
                surf(dgaussh);
                handles.J=dgaussh;
                guidata(hObject, handles)

                axes(handles.axes4);
                imshow(bgaussh);
                handles.J=bgaussh;
                guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of gaussian_filter
% --- Executes on button press in butter_filter.
        else 
            m=220;   % Size is set
            n=220;
            d0=handles.o;
            a=im2double(rgb2gray(imread('lena.png')));
            N=1;                        % Set the filter order
            % Determine the centre of the image
            p= round(m/2);
            q= round(n/2);
                %Define the kernel
            H= zeros(m,n);
            for i = 1:m
                for j = 1:n
                    d = (i-p).^2+(j-q).^2;
                    if d~=0
                    H(i,j) = 1/(1+((d0*d0/d).^(2*N)));
                    end
                end
            end
            %Input image in frequency domain
            A_f = fftshift(fft2(a));
            %Apply the filter
            B = A_f.*H;
            C = abs(ifft2(B));
            %Display
            axes(handles.axes2);
            imshow(H);
            handles.J=H;
            guidata(hObject, handles);

            axes(handles.axes3);
            surf(H);
            handles.J=H;
            guidata(hObject, handles);

            axes(handles.axes4);
            imshow(C);
            handles.J=C;
            guidata(hObject, handles);

        end
 
    end


% --- Executes during object creation, after setting all properties.
function high_pass_CreateFcn(hObject, eventdata, handles)
set(hObject,'Value',1); % if you want it to be seelcted initialy
% hObject    handle to low_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
function show_image_Callback(hObject, eventdata, handles)
% hObject    handle to show_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=im2double(rgb2gray(imread('lena.png')));

% Convert to grayscale incase it is color
imshow(a,'Parent',handles.axes1) ;




% Hint: get(hObject,'Value') returns toggle state of high_pass

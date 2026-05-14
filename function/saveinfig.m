function saveinfig(f,imgTitle,tag,subdir)
if nargin<4 ||isempty(subdir)
    st = dbstack;
%     st(1).name % The function's name
     % The function caller's name (parent)
    subdir=st(2).name;
end
if ~iscell(f)
    f={f};
end
if ischar(imgTitle)
    imgTitle={imgTitle};
end

for i = 1:length(f)
    img=f{i};
    imageName = ['./fig/',subdir,filesep,tag,filesep, imgTitle{i}, '.fig'];
    test_folder(imageName);
    saveas(img, imageName);
    imageName = ['./fig/',subdir,filesep,tag,filesep, imgTitle{i}, '.jpg'];
    test_folder(imageName);
    saveas(img, imageName);
end
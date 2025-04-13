function VisualizeReg_MDACCbcdata_newSet(subname,matpath,i)

for v = 1:3
   filename = [matpath subname, '_v',num2str(v),'_Reg.mat'];
   V = load(filename);
   eval(['V',num2str(v),' = V;'])
   if i == 1
    vuOnePaneViewer(V.anatomical)
    vuOnePaneViewer(V.roi)
   end
end

vuOnePaneViewer(V1.anatomical+V2.anatomical+V3.anatomical)
vuOnePaneViewer(V1.roi+10*V2.roi+20*V3.roi)

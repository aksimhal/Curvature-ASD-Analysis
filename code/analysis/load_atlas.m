function atlas_labels = load_atlas()

fileID = fopen('atlas_key.csv','r');
C = textscan(fileID,'%s %s');
atlas_labels = C{2};

for n =1:length(atlas_labels)
    atlas_labels{n} = atlas_labels{n}(1:end-1);
end

end
    
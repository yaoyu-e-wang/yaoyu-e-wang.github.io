#! /bin/bash

# assuming Hugo has been install locally
site=$1

echo "creating hugo website: $1"

# create a new site
hugo new site $site

# Add a Theme

cd $site
git init

# Default adding ananke theme
# git submodule add https://github.com/budparr/gohugo-theme-ananke.git themes/ananke

# Adding ramium theme
git submodule add https://github.com/rafed/ramium.git themes/ramium
cp ../templates/ramium/config.toml config.toml

# Adding Even theme and copy config.toml to root 
# git clone https://github.com/olOwOlo/hugo-theme-even themes/even
#cp ../templates/even/config.toml config.toml

# use customized color theme at themes/even/assets/sass/_variables.scss
cp ../templates/ramium/_variables.scss themes/ramium/assets/sass/_variables.scss

#git submodule update --remote --merge

# add theme into config file
echo 'theme = "ramium"' >> config.toml

# Add test posts for anake and ramium use 'posts'
hugo new posts/my-test-post.md
echo 'This is the test post' > content/posts/my-first-post.md
# hugo new posts/my-second-post.md
# cat content/posts/my-second-post.md ../example_md/README.md > content/posts/my-second-post.md
# cp -r ../example_md/*.md content/posts/

# Add posts for even, use 'post'
# hugo new post/my-first-post.md
# echo 'This is the first post' >> content/post/my-first-post.md

# create content post directory
# mkdir content/post/
# cp -r ../notes_markdown/*.md content/post/
cp -r ../notes_markdown/*.md content/post/
cp -r ../notes_markdown/images/ static/images

# start Hugo server with draft posts enabled

hugo server -D

# build static pages:

hugo -D

# Output will be in ./public/ directory by default 
# (-d/--destination flag to change it, or set publishdir in the config file).







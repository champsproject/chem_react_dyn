N=$1
rm -rf _build _site
jupyter-book build ./
bundle exec jekyll serve --watch --port $N

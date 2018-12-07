#
echo did you refresh the pages by running?
echo 'julia --color=yes make.jl && cp -a build/* site/'
echo 'firefox site/index.html'
echo "press enter to continue. Ctrl-C to stop"
read -n 1
echo 'Publishing'
export DOCSDIR=`pwd`
export newdir=/tmp/neppack$$
mkdir $newdir
cd $newdir
git clone -b "gh-pages" git@github.com:nep-pack/NonlinearEigenproblems.jl.git
cd NonlinearEigenproblems.jl/
cp -r $DOCSDIR/site/* .
git add *;  git commit . -m "refresh docs"; git push
cd $DOCSDIR
echo "rm -Rf $newdir"

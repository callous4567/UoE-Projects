#**This directory just has some old deprecated/useless codes.** 
They may be useful as reference material, hence they've been retained. 
A few are fully-functioning implementations of some papers. 
By no means is the code here entirely worth discarding, but they should *not* be taken as being "fully functional." 
You should also note that *windows_testing* codes are ones that have not been fully functional at any point as part
of the work pipeline. 


###*windows_confimap* 
For example, is part of an implementation of Charles King III use of Great Circle Cell counting. 
While fun and interesting, and capable at reproducing the maps in that paper, this code proved futile for the date 
I had to play with. Mind you I did try it on the original dataset and it produced reasonable results. The reason it 
it failed here was that as n-phi got bigger and bigger, you get to the point where an entire great-circle of poles 
becomes significant for a given region of space. 

Envision a spatial clump and draw a GCC with n-phi=1 (i.e. the full GCC.) In this instance, your GCC count will yield
one reasonable pole- all the others aren't useful. If instead you defined n-phi as some arbitrarily large value, 
you could consider rotating that original GCC plane about the vector from the centre of the galcentric frame
to the spatial clump- no matter the rotation, all the GCCs have the same confidence, since they're all using the same 
area for the same clump that overlaps that area- i.e. the same "significant" number density. In any case, this 
manifests as long streaks (great circle streaks) in the final "confidence map" which are totally worthless to the final
analysis.

Now, you might ask why that was a problem- well it's a real big one, since our data is pretty lush (~30,000 stars, 
if you consider all the sets as one as planned.) Consider trying to get the binomial CDF for very large numbers- 
it. does. not. work. You end up with infinities when you decide to do the map, log(1/p) ends up nearing infinity, and
the result is white space in your confidence map- larger cell areas yields this effect- yet we need larger cell areas 
to actually resolve the finer detail needed to get the poles in the data, which ends up introducing the null-result of
great circles of GCCs on the map as mentioned earlier. It's a catch-22... a balancing act... and one that is futile
for our data. Thus it was abandoned. In any case, the code is fully-functional and useful, so if you want it, feel free.

You could *improve and make it more robust* to this problem by running the routine as normal, but in the stage of 
grabbing the "most confident" phi-subcell for the mapping from each GCC pole, select it, rotate to center on that cell,
then try re-running. This would let you average out any spontaneous coincidences in just doing [0, 360/nphi,...]. 
**Note that setting nphi to 1 just returns you standard Johnston 1996 Great Circle Cell Counting**. 


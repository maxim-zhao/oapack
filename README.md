# Optimal aPACK compressor

*aPACK* is data compression format developed by [Jørgen Ibsen](http://ibsensoftware.com/products_aPLib.html). Here is some specific compressor implementation in C language.

The idea behind this particular compressor is achieving maximum compression possible for the given compression scheme. It does exhaustive search, trying all possible states using dynamic programming. For this reason, the compressor is very slow and has the following **limitations**:

*  input file size may not exceed 65536 bytes;
*  file must not contain long repetitive runs; processing these is extremely slow;
*  severe amounts of memory, up to 9\**N*<sup>2</sup> bytes required; for this reason 64-bit mode is recommended.

The implementation is pretty straight and does not involve complex data structures. The resulting algorithm complexity is *O*(*n*<sup>3</sup>), which seems insane, but in most practical cases it finishes in reasonable time span.

This compressor produces raw compressed stream, no headers, no decompressor included.

Supports back-to-front compression (*-b* option).

### Other aPACK compressors

* [aplib_pack2 by r57shell](http://gendev.spritesmind.net/forum/viewtopic.php?f=7&t=703&&start=45#p32548)
* [apc12spke](https://www.cpcwiki.eu/forum/programming/quick-update-on-the-state-of-the-art-compression-using-aplib/msg177112/#msg177112)
* [apultra](https://github.com/emmanuel-marty/apultra)

function ds = downsample(s, factor)
  % downsample a signal by factor. the original and downsampled signals
  % have equal integrals.
  %
  %        $Id$
  % to downsample conserving the amplitude use: 
  %     ds = downsample(s, factor) / factor;
  
  if size(s,1)==1 %if only one row, we've been passed a row vector
     transpose_vector = 1;
     s = s';
  else
     transpose_vector = 0;
  end
  
  a = cumsum(s, 1);
  ds = diff([zeros(1, size(s,2));a(factor:factor:end,:)]);

  if transpose_vector
     ds = ds';
  end
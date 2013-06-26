      program createheader

      open(unit=2,file='header.bin',form='unformatted')
      write(2) 1
      write(2) 1.0d0
      write(2) .true.
      close(2)
      end

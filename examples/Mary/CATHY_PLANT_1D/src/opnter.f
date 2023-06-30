      subroutine opnter(luname)
      character term*20,random*20,comand*80
      comand='tty > /tmp/open'
      call system(comand)
      open(99,file='/tmp/open')
      read (99,'(A)') term
      close(99,status='delete')
      open(luname,file=term)
      return
      end

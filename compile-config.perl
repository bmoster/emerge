######################################################################################
# Emerge code - compile-config.perl                                                  #
# Adapted from the GADGET code developed by Volker Springel                          #
######################################################################################
#                                                                                    #
# This file processes the configurations options in Config.sh, producing two files:  #
#   codeoptions.h        to be included in each source file (via allvars.h)          #
#   compile_info.c       code to be compiled in, which will print the configuration  #
#                                                                                    #
######################################################################################

if( @ARGV != 2)
{
  print "usage: perl compile-config.perl <Config.sh> <build dir>\n";
  exit;
}

open(FILE, @ARGV[0]);
$path = @ARGV[1];

open(OUTFILE, ">${path}/codeoptions.h");
open(COUTF,   ">${path}/compile_info.c");

print COUTF "#include <stdio.h>\n";
print COUTF "#include \"../src/allvars.h\"\n\n";
print COUTF "void output_code_options(void)\n\{\n";

$count = 0;

while($line=<FILE>)
{
  chop $line;

  @fields = split ' ' , $line;

  if(substr($fields[0], 0, 1) ne "#")
  {
    if(length($fields[0]) > 0)
    {
      @subfields = split '=', $fields[0];

      print OUTFILE "#define $subfields[0] $subfields[1]\n";
      print COUTF   "    printf(\"%c%c%c $fields[0]\\n\",SYMBOL,SYMBOL,SYMBOL);\n";
      $count = $count + 1;
    }
  }
}

if ($count<1)
{
  print COUTF "    printf(\"%c%c%c NONE\\n\",SYMBOL,SYMBOL,SYMBOL);\n";
}

print COUTF "\}\n";

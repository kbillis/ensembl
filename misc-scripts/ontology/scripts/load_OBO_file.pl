#!/usr/bin/perl -w

# A simple OBO file reader/loader

use strict;
use warnings;

use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );
use IO::File;

#-----------------------------------------------------------------------

sub usage {
  print("Usage:\n");
  printf( "\t%s\t-h dbhost [-P dbport] \\\n"
      . "\t%s\t-u dbuser [-p dbpass] \\\n"
      . "\t%2\$s\t-d dbname [-t] \\\n"
      . "\t%2\$s\t-f file\n",
    $0, ' ' x length($0) );
  print("\n");
  printf( "\t%s\t-?\n", $0 );
  print("\n");
  print("Arguments:\n");
  print("\t-h/--host dbhost\tDatabase server host name\n");
  print("\t-P/--port dbport\tDatabase server port (optional)\n");
  print("\t-u/--user dbuser\tDatabase user name\n");
  print("\t-p/--pass dbpass\tUser password (optional)\n");
  print("\t-d/--name dbname\tDatabase name\n");
  print("\t-t/--truncate\t\tTruncate (empty) each table\n");
  print("\t\t\t\tbefore writing (optional)\n");
  print("\t-f/--file file\t\tThe OBO file to parse\n");
  print("\t-?/--help\t\tDisplays this information\n");
}

#-----------------------------------------------------------------------

sub write_ontology {
  my ( $dbh, $truncate, $namespaces ) = @_;

  print("Writing to 'ontology' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE ontology") }

  my $statement = "INSERT INTO ontology (name, namespace) VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $namespace ( sort( keys( %{$namespaces} ) ) ) {
    my $ontology = $namespaces->{$namespace};

    $sth->bind_param( 1, $ontology,  SQL_VARCHAR );
    $sth->bind_param( 2, $namespace, SQL_VARCHAR );

    $sth->execute();

    $namespaces->{$namespace} = {
      'id' =>
        $dbh->last_insert_id( undef, undef, 'ontology', 'ontology_id' ),
      'name' => $ontology
    };

    ++$count;
  }

  $dbh->do("OPTIMIZE TABLE ontology");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_ontology

#-----------------------------------------------------------------------

sub write_subset {
  my ( $dbh, $truncate, $subsets ) = @_;

  print("Writing to 'subset' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE subset") }

  $dbh->do("LOCK TABLE subset WRITE");

  my $statement =
    "INSERT INTO subset " . "(name, definition) " . "VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $subset_name ( sort( keys( %{$subsets} ) ) ) {
    my $subset = $subsets->{$subset_name};

    $sth->bind_param( 1, $subset->{'name'},       SQL_VARCHAR );
    $sth->bind_param( 2, $subset->{'definition'}, SQL_VARCHAR );

    $sth->execute();

    $subset->{'id'} =
      $dbh->last_insert_id( undef, undef, 'subset', 'subset_id' );
    ++$count;
  }

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_subset

#-----------------------------------------------------------------------

sub write_term {
  my ( $dbh, $truncate, $terms, $subsets, $namespaces ) = @_;

  print("Writing to 'term' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE term") }

  $dbh->do("LOCK TABLE term WRITE");

  my $statement =
      "INSERT INTO term "
    . "(ontology_id, subsets, accession, name, definition) "
    . "VALUES (?,?,?,?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $accession ( sort( keys( %{$terms} ) ) ) {
    my $term = $terms->{$accession};

    my $term_subsets;

    if ( @{ $term->{'subsets'} } ) {
      $term_subsets = join( ',',
        map { $subsets->{$_}{'name'} } @{ $term->{'subsets'} } );
    }

    $sth->bind_param( 1, $namespaces->{ $term->{'namespace'} }{'id'},
      SQL_INTEGER );
    $sth->bind_param( 2, $term_subsets,         SQL_VARCHAR );
    $sth->bind_param( 3, $accession,            SQL_VARCHAR );
    $sth->bind_param( 4, $term->{'name'},       SQL_VARCHAR );
    $sth->bind_param( 5, $term->{'definition'}, SQL_VARCHAR );

    $sth->execute();

    $term->{'id'} =
      $dbh->last_insert_id( undef, undef, 'term', 'term_id' );
    ++$count;
  }

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_term

#-----------------------------------------------------------------------

sub write_relation_type {
  my ( $dbh, $truncate, $relation_types ) = @_;

  print("Writing to 'relation_type' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE relation_type") }

  my $select_stmt =
    "SELECT relation_type_id FROM relation_type WHERE name = ?";
  my $select_sth = $dbh->prepare($select_stmt);

  my $insert_stmt = "INSERT INTO relation_type (name) VALUES (?)";
  my $insert_sth  = $dbh->prepare($insert_stmt);

  my $count = 0;
  foreach my $relation_type ( sort( keys( %{$relation_types} ) ) ) {
    $select_sth->bind_param( 1, $relation_type, SQL_VARCHAR );
    $select_sth->execute();

    my $id;
    my $found = 0;

    $select_sth->bind_columns( \$id );

    while ( $select_sth->fetch() ) {
      $relation_types->{$relation_type} = { 'id' => $id };
      $found = 1;
    }

    if ( !$found ) {
      $insert_sth->bind_param( 1, $relation_type, SQL_VARCHAR );
      $insert_sth->execute();
      $relation_types->{$relation_type} = {
        'id' => $dbh->last_insert_id(
          undef, undef, 'relation_type', 'relation_type_id'
        ) };
      ++$count;
    }
  }

  $dbh->do("OPTIMIZE TABLE relation_type");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_relation_type

#-----------------------------------------------------------------------

sub write_relation {
  my ( $dbh, $truncate, $terms, $relation_types ) = @_;

  print("Writing to 'relation' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE relation") }

  $dbh->do("LOCK TABLE relation WRITE");

  my $statement =
      "INSERT INTO relation "
    . "(child_term_id, parent_term_id, relation_type_id) "
    . "VALUES (?,?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $child_term ( sort { $a->{'id'} <=> $b->{'id'} }
    values( %{$terms} ) )
  {
    foreach my $relation_type (
      sort( keys( %{ $child_term->{'parents'} } ) ) )
    {
      foreach my $parent_acc (
        sort( @{ $child_term->{'parents'}{$relation_type} } ) )
      {
        $sth->bind_param( 1, $child_term->{'id'},         SQL_INTEGER );
        $sth->bind_param( 2, $terms->{$parent_acc}{'id'}, SQL_INTEGER );
        $sth->bind_param( 3, $relation_types->{$relation_type}{'id'},
          SQL_INTEGER );

        $sth->execute();

        ++$count;
      }
    }
  }

  $dbh->do("OPTIMIZE TABLE relation");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_relation

#-----------------------------------------------------------------------

my ( $dbhost, $dbport );
my ( $dbuser, $dbpass );
my ( $dbname, $truncate, $obo_file_name );

$dbport   = '3306';
$truncate = 0;

if (
  !GetOptions(
    'dbhost|host|h=s' => \$dbhost,
    'dbport|port|P=i' => \$dbport,
    'dbuser|user|u=s' => \$dbuser,
    'dbpass|pass|p=s' => \$dbpass,
    'dbname|name|d=s' => \$dbname,
    'truncate|t'      => \$truncate,
    'file|f=s'        => \$obo_file_name,
    'help|?'          => sub { usage(); exit } )
  || !defined($dbhost)
  || !defined($dbuser)
  || !defined($dbname)
  || !defined($obo_file_name) )
{
  usage();
  exit;
}

my $dsn = sprintf( 'dbi:mysql:database=%s;host=%s;port=%s',
  $dbname, $dbhost, $dbport );

my $dbh =
  DBI->connect( $dsn, $dbuser, $dbpass,
  { 'RaiseError' => 1, 'PrintError' => 2 } );

my $statement =
  "SELECT meta_value FROM meta WHERE meta_key = 'OBO_file_date'";

my $stored_obo_file_date = $dbh->selectall_arrayref($statement)->[0][0];
my $obo_file_date;

my $obo = IO::File->new( $obo_file_name, 'r' ) or die;

my ( $state, $accession, $name, $namespace, $definition );
my ( %parents, @subsets );

my ( %terms, %namespaces, %relation_types, %subsets );

printf( "Reading OBO file '%s'...\n", $obo_file_name );

my $default_namespace;

while ( defined( my $line = $obo->getline() ) ) {
  chomp($line);

  if ( !defined($state) ) {
    if ( $line =~ /^\[(\w+)\]$/ ) { $state = $1; next }

    if ( $line =~ /^([\w-]+): (.+)$/ ) {
      my $type = $1;
      my $data = $2;

      if ( $type eq 'date' ) {
        $obo_file_date = sprintf( "%s/%s", $obo_file_name, $data );

        if ( defined($stored_obo_file_date) ) {
          if ( $stored_obo_file_date eq $obo_file_date
            && !$truncate )
          {
            print("This OBO file has already been processed.\n");
            $obo->close();
            $dbh->disconnect();
            exit;
          } elsif ( index( $stored_obo_file_date, $obo_file_name ) != -1
            && !$truncate )
          {
            print <<EOT;
==> Trying to load a newer (?) OBO file that has already been loaded.
==> Please clean the database manually of data associated with this
==> file and try again... or use the -t (truncate) switch to empty the
==> tables completely (unless you want to preserve some of the data,
==> obviously).
EOT
            $obo->close();
            $dbh->disconnect();
            exit;
          }
        }

      } elsif ( $type eq 'default-namespace' ) {
        $default_namespace = $data;
      } elsif ( $type eq 'subsetdef' ) {
        my ( $subset_name, $subset_def ) =
          ( $data =~ /^(\w+)\s+"([^"]+)"/ );
        $subsets{$subset_name}{'name'}       = $subset_name;
        $subsets{$subset_name}{'definition'} = $subset_def;
      }
    } ## end if ( $line =~ /^([\w-]+): (.+)$/)

    next;
  } ## end if ( !defined($state) )

  if ( $state eq 'Term' ) {
    if ( $line eq '' ) {
      $namespace ||= $default_namespace;
      ( $namespaces{$namespace} ) = $accession =~ /^([^:]+):/;

      $terms{$accession} = {
        'namespace'  => $namespace,
        'name'       => $name,
        'definition' => $definition,
        'subsets'    => [@subsets] };

      foreach my $relation_type ( keys(%parents) ) {
        if ( !exists( $relation_types{$relation_type} ) ) {
          $relation_types{$relation_type} = 1;
        }
        $terms{$accession}{'parents'} = {%parents};
      }

      $state = 'clear';
    } elsif ( $line =~ /^(\w+): (.+)$/ ) {
      my $type = $1;
      my $data = $2;

      if    ( $type eq 'id' )        { $accession = $data }
      elsif ( $type eq 'name' )      { $name      = $data }
      elsif ( $type eq 'namespace' ) { $namespace = $data }
      elsif ( $type eq 'def' ) {
        ($definition) = $data =~ /"([^"]+)"/;
      } elsif ( $type eq 'is_a' ) {
        my ($parent_acc) = $data =~ /(\S+)/;
        push( @{ $parents{'is_a'} }, $parent_acc );
      } elsif ( $type eq 'relationship' ) {
        my ( $relation_type, $parent_acc ) = $data =~ /^(\w+) (\S+)/;
        push( @{ $parents{$relation_type} }, $parent_acc );
      } elsif ( $type eq 'is_obsolete' ) {
        if ( $data eq 'true' ) { $state = 'clear' }
      } elsif ( $type eq 'subset' ) {
        push( @subsets, $data );
      }

    }
  } ## end if ( $state eq 'Term' )

  if ( $state eq 'clear' ) {
    undef($accession);
    undef($name);
    undef($namespace);
    undef($definition);
    %parents = ();
    @subsets = ();

    undef($state);
  }

} ## end while ( defined( my $line...))

$obo->close();

print("Finished reading OBO file, now writing to database...\n");

write_ontology( $dbh, $truncate, \%namespaces );
write_subset( $dbh, $truncate, \%subsets );
write_term( $dbh, $truncate, \%terms, \%subsets, \%namespaces );
write_relation_type( $dbh, $truncate, \%relation_types );
write_relation( $dbh, $truncate, \%terms, \%relation_types );

print("Updating meta table...\n");

if ($truncate) {
  $dbh->do("TRUNCATE TABLE meta");
}

my $sth =
  $dbh->prepare( "DELETE FROM meta "
    . "WHERE meta_key = 'OBO_file_date' "
    . "AND meta_value LIKE ?" );
$sth->bind_param( 1, sprintf( "%s/%%", $obo_file_name ), SQL_VARCHAR );
$sth->execute();
$sth->finish();

$sth = $dbh->prepare( "INSERT INTO meta (meta_key, meta_value)"
    . "VALUES ('OBO_file_date', ?)" );
$sth->bind_param( 1, $obo_file_date, SQL_VARCHAR );
$sth->execute();
$sth->finish();

my $obo_load_date =
  sprintf( "%s/%s", $obo_file_name, scalar( localtime() ) );

$sth =
  $dbh->prepare( "DELETE FROM meta "
    . "WHERE meta_key = 'OBO_load_date' "
    . "AND meta_value LIKE ?" );
$sth->bind_param( 1, sprintf( "%s/%%", $obo_file_name ), SQL_VARCHAR );
$sth->execute();
$sth->finish();

$sth = $dbh->prepare( "INSERT INTO meta (meta_key, meta_value)"
    . "VALUES ('OBO_load_date', ?)" );
$sth->bind_param( 1, $obo_load_date, SQL_VARCHAR );
$sth->execute();
$sth->finish();

$dbh->disconnect();

print("Done.\n");

# $Id$

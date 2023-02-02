// write a code to reverse the array?
void reverse(char word[])
{
    int len=strlen(word);
    char temp;
    for (int i=0;i<len/2;i++)
    {
            temp=word[i];
            word[i]=word[len-i-1];
            word[len-i-1]=temp;
    }
}


void reverse(char word[])
{
    int len=strlen(word);
    for (int i=0;i<len/2;i++)
    {
        word[i]^=word[len-i-1];
        word[len-i-1]^=word[i];
        word[i]^=word[len-i-1];
    }
}




//Source: https://stackoverflow.com/questions/1128985



r = rand() % 11;


r = rand() % (n + 1);


r = rand() % (n + 1) + k;




//Source: https://stackoverflow.com/questions/4919303



